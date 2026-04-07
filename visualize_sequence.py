#!/usr/bin/env python3
"""
Interactive DNA sequence viewer with SnapGene-like double-stranded display.

Features:
  * Forward strand on top, complement strand beneath (aligned per-base).
  * Interactive search bar with live highlighting (both strands).
  * Navigate between matches with Prev/Next buttons (or Enter/Shift+Enter).
  * Sequence selector dropdown for multi-sequence FASTA files.
  * Virtual scrolling — renders only visible blocks for instant performance.
  * Optional BLAST pre-rendering from CLI (green highlights persist alongside search).
  * Live BLAST from the GUI when running in --serve mode.
  * "Go to position" field to jump anywhere in the sequence.

Usage:
  python visualize_sequence.py <fasta> [options]

Examples:
  # Static HTML (open in browser / Live Preview)
  python visualize_sequence.py references/6991.fasta -o 6991.html

  # Interactive server with live BLAST from the GUI
  python visualize_sequence.py references/6991.fasta --serve

  # Pre-baked BLAST hits
  python visualize_sequence.py references/6991.fasta --blast query.fa -o 6991.html
"""

import argparse
import html
import json
import os
import subprocess
import sys
import tempfile
import webbrowser
from http.server import HTTPServer, BaseHTTPRequestHandler
from pathlib import Path

COMPLEMENT = str.maketrans("ACGTNRYSWKMBDHVacgtnryswkmbdhv",
                           "TGCANYRSWMKVHDBtgcanyrswmkvhdb")


def parse_fasta(path):
    records = []
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(buf).upper()))
                name = line[1:].split()[0]
                buf = []
            elif line:
                buf.append(line)
        if name is not None:
            records.append((name, "".join(buf).upper()))
    if not records:
        sys.exit(f"No sequences found in {path}")
    return records


def complement_str(seq):
    return seq.translate(COMPLEMENT)


def run_blast(subject_fasta, query_arg, evalue=10.0, word_size=7):
    if not any(os.access(os.path.join(p, "blastn"), os.X_OK)
               for p in os.environ.get("PATH", "").split(os.pathsep)):
        raise RuntimeError("blastn not found in PATH.")

    cleanup = []
    try:
        if os.path.isfile(query_arg):
            query_path = query_arg
        else:
            tf = tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False)
            tf.write(f">query\n{query_arg}\n")
            tf.close()
            query_path = tf.name
            cleanup.append(query_path)

        cmd = [
            "blastn", "-subject", subject_fasta, "-query", query_path,
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                       "qstart qend sstart send evalue bitscore sstrand",
            "-evalue", str(evalue), "-word_size", str(word_size), "-dust", "no",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    finally:
        for p in cleanup:
            try:
                os.unlink(p)
            except OSError:
                pass

    hits = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        f = line.split("\t")
        sstart, send = int(f[8]), int(f[9])
        strand = "+" if f[12] == "plus" else "-"
        if sstart > send:
            sstart, send = send, sstart
        hits.append({
            "id": f[1], "pi": round(float(f[2]), 1), "l": int(f[3]),
            "s": sstart - 1, "e": send,
            "ev": float(f[10]), "bs": round(float(f[11]), 1), "st": strand,
        })
    return hits


# ── HTML / CSS / JS ──────────────────────────────────────────────────

CSS = r"""
* { box-sizing: border-box; margin: 0; padding: 0; }
html, body { height: 100%; overflow: hidden; }
body {
    background: #1e1e1e; color: #d4d4d4;
    font-family: 'Courier New', Consolas, monospace;
    display: flex; flex-direction: column;
}

/* ── toolbar ── */
#toolbar {
    flex: none; z-index: 100;
    background: #252526; border-bottom: 1px solid #3c3c3c;
    padding: 8px 16px; display: flex; align-items: center; gap: 12px;
    flex-wrap: wrap;
}
#toolbar label { font-size: 12px; color: #aaa; }
#toolbar select, #toolbar input[type="text"] {
    background: #3c3c3c; color: #e0e0e0; border: 1px solid #555;
    padding: 4px 8px; font-size: 13px; font-family: inherit;
    border-radius: 3px;
}
#toolbar select { max-width: 220px; }
#searchInput { width: 260px; }
#gotoInput { width: 100px; }
#toolbar button {
    background: #0e639c; color: #fff; border: none; padding: 4px 12px;
    border-radius: 3px; cursor: pointer; font-size: 12px;
}
#toolbar button:hover { background: #1177bb; }
#toolbar button:disabled { opacity: 0.4; cursor: default; }
#matchInfo { font-size: 12px; color: #aaa; min-width: 90px; }
#toolbar .sep { width: 1px; height: 24px; background: #555; }
#rcLabel { font-size: 12px; color: #aaa; display: flex; align-items: center; gap: 4px; cursor: pointer; }
#rcLabel input { accent-color: #0e639c; }
#posIndicator { font-size: 12px; color: #666; margin-left: auto; }

/* ── info bar (legend, blast summary) ── */
#infobar { flex: none; padding: 8px 24px; background: #1e1e1e; }
.legend { font-size: 12px; color: #bbb; margin-bottom: 4px; }
.legend .swatch { display: inline-block; padding: 1px 6px; margin-right: 4px; border-radius: 2px; font-weight: bold; }
.blast-summary {
    background: #2a2a2a; padding: 8px 12px; border-left: 3px solid #569cd6;
    margin-top: 6px; font-size: 12px; max-height: 150px; overflow-y: auto;
}
.blast-summary table { border-collapse: collapse; font-size: 12px; }
.blast-summary td, .blast-summary th { padding: 2px 10px 2px 0; text-align: left; }
.blast-summary th { color: #9cdcfe; font-weight: normal; }

/* ── viewport (virtual scroll) ── */
#viewport {
    flex: 1; overflow-y: auto; overflow-x: hidden;
    position: relative;
}
#spacer { position: relative; width: 100%; user-select: none; cursor: text; }

/* ── sequence blocks ── */
.block {
    position: absolute; left: 0; right: 0;
    padding: 0 24px;
}
.c { /* base char wrapper */ }
.sel { background: #264f78 !important; color: #fff !important; font-weight: normal !important; }
.line { white-space: pre; font-size: 13px; line-height: 1.4; letter-spacing: 1px; }
.ruler { color: #555; font-size: 11px; letter-spacing: 1px; line-height: 1.3; }
.fwd { color: #e0e0e0; padding-bottom: 1px; }
.ticks { color: #444; font-size: 11px; letter-spacing: 1px; line-height: 0.8; }
.rev { color: #7a7a7a; padding-top: 1px; }
.num { color: #569cd6; display: inline-block; min-width: 56px;
       text-align: right; margin-right: 8px; letter-spacing: 0; }

.hl { font-weight: bold; border-radius: 2px; padding: 0 1px; }
.hl-search { background: #ffd54a; color: #111; }
.hl-blast  { background: #4caf50; color: #111; }
.hl-active { background: #ff6d00; color: #111; }

/* ── blast panel ── */
#blastInput { width: 260px; }
#blastInfo { font-size: 12px; color: #aaa; }
#blastSettingsBtn { background: transparent; color: #aaa; border: 1px solid #555;
    font-size: 14px; padding: 2px 7px; cursor: pointer; }
#blastSettingsBtn:hover { color: #e0e0e0; border-color: #888; }
#blastSettingsBtn.active { color: #4fc3f7; border-color: #4fc3f7; }

/* ── blast settings panel ── */
#blastSettings {
    display: none; flex: none; background: #2a2a2a; border-bottom: 1px solid #3c3c3c;
    padding: 8px 24px; gap: 16px; align-items: center; font-size: 12px; color: #aaa;
    flex-wrap: wrap;
}
#blastSettings.open { display: flex; }
#blastSettings label { display: flex; align-items: center; gap: 4px; }
#blastSettings input[type="number"] {
    background: #3c3c3c; color: #e0e0e0; border: 1px solid #555;
    padding: 3px 6px; width: 70px; font-size: 12px; font-family: inherit;
    border-radius: 3px;
}
#blastSettings .filter-title { color: #9cdcfe; margin-right: 4px; }
#filteredCount { color: #888; font-style: italic; }

/* ── position tooltip ── */
#posTip {
    position: fixed; z-index: 200; pointer-events: none;
    background: #1e1e1e; color: #ccc; border: 1px solid #555;
    padding: 3px 8px; border-radius: 3px; font-size: 12px;
    font-family: 'Courier New', Consolas, monospace;
    white-space: nowrap; display: none;
}
.blast-running { opacity: 0.6; pointer-events: none; }
"""

JS = r"""
const W = LINE_WIDTH;
const BUFFER = 30;

let curIdx = 0;
let matches = [];
let activeMatch = -1;
let lastQuery = '';
let searchMarks = null;   // Uint8Array over current seq
let blastCls = null;      // Uint8Array over current seq
let numBlocks = 0;
let BLOCK_H = 0;
let renderedBlocks = {};  // blockIdx -> DOM element

const seqSel    = document.getElementById('seqSelect');
const searchIn  = document.getElementById('searchInput');
const matchInfo = document.getElementById('matchInfo');
const prevBtn   = document.getElementById('prevBtn');
const nextBtn   = document.getElementById('nextBtn');
const gotoIn    = document.getElementById('gotoInput');
const rcCheck   = document.getElementById('rcCheck');
const viewport  = document.getElementById('viewport');
const spacer    = document.getElementById('spacer');
const posInd    = document.getElementById('posIndicator');

// ── complement ──
const COMP = {};
'ACGTNRYSWKMBDHVacgtnryswkmbdhv'.split('').forEach((c,i) => {
    COMP[c] = 'TGCANYRSWMKVHDBtgcanyrswmkvhdb'[i];
});
function rc(s) {
    let r = '';
    for (let i = s.length - 1; i >= 0; i--) r += (COMP[s[i]] || s[i]);
    return r;
}

// ── measure block height ──
function measureBlockH() {
    const div = document.createElement('div');
    div.className = 'block';
    div.style.position = 'static';
    div.style.visibility = 'hidden';
    div.innerHTML =
        '<div class="line ruler"><span class="num">&nbsp;</span>' + '\u00b7'.repeat(10) + '</div>' +
        '<div class="line fwd"><span class="num">&nbsp;</span>ACGTACGTAC</div>' +
        '<div class="line ticks"><span class="num">&nbsp;</span>' + '\u2500'.repeat(10) + '</div>' +
        '<div class="line rev"><span class="num">&nbsp;</span>TGCATGCATG</div>';
    document.body.appendChild(div);
    BLOCK_H = Math.ceil(div.getBoundingClientRect().height) + 8; // +8 for gap
    document.body.removeChild(div);
}

// ── blast marks + navigation ──
let blastHitsCur = [];    // hits on current seq, sorted by position
let activeBlastHit = -1;

function buildBlastCls(rec) {
    const n = rec.seq.length;
    blastCls = new Uint8Array(n);
    const entry = BLAST_MARKS.find(b => b.id === rec.id);
    if (entry) {
        for (const m of entry.ranges) {
            for (let i = m.s; i < m.e && i < n; i++) blastCls[i] = 1;
        }
    }
    // mark active blast hit
    if (activeBlastHit >= 0 && activeBlastHit < blastHitsCur.length) {
        const h = blastHitsCur[activeBlastHit];
        for (let i = h.s; i < h.e && i < n; i++) blastCls[i] = 2;
    }
}

function buildBlastHitsCur() {
    const rec = SEQS[curIdx];
    const entry = BLAST_MARKS.find(b => b.id === rec.id);
    blastHitsCur = entry ? entry.ranges.slice().sort((a, b) => a.s - b.s) : [];
    activeBlastHit = blastHitsCur.length > 0 ? 0 : -1;
}

// ── render a single block ──
function createBlockEl(rec, bi) {
    const start = bi * W;
    const end = Math.min(start + W, rec.seq.length);
    const seq = rec.seq;
    const comp = rec.comp;

    const sMin = Math.min(selStart, selEnd);
    const sMax = Math.max(selStart, selEnd);
    const hasSel = selStart >= 0;

    let ruler = '', fwd = '', rev = '', ticks = '';
    for (let i = start; i < end; i++) {
        const col = i - start + 1;
        ruler += (col % 10 === 0) ? '|' : '\u00b7';
        ticks += (col % 10 === 0) ? '|' : '\u2500';

        const sm = searchMarks ? searchMarks[i] : 0;
        const bc = blastCls[i];
        let cls = 'c';
        if (sm === 2) cls = 'c hl hl-active';
        else if (sm === 1) cls = 'c hl hl-search';
        else if (bc === 2) cls = 'c hl hl-active';
        else if (bc === 1) cls = 'c hl hl-blast';

        if (hasSel && i >= sMin && i <= sMax) cls += ' sel';

        fwd += '<span data-p="' + i + '" class="' + cls + '">' + seq[i] + '</span>';
        rev += '<span data-p="' + i + '" class="' + cls + '">' + comp[i] + '</span>';
    }

    const div = document.createElement('div');
    div.className = 'block';
    div.style.top = (bi * BLOCK_H) + 'px';
    div.dataset.bi = bi;

    const lbl = String(start + 1).padStart(7, '\u00a0');
    div.innerHTML =
        '<div class="line ruler"><span class="num">\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0</span>' + ruler + '</div>' +
        '<div class="line fwd"><span class="num">' + lbl + '</span>' + fwd + '</div>' +
        '<div class="line ticks"><span class="num">\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0</span>' + ticks + '</div>' +
        '<div class="line rev"><span class="num">\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0</span>' + rev + '<span class="num"> ' + end + '</span></div>';
    return div;
}

// ── virtual scroll: render/remove blocks ──
function updateVisible() {
    const scrollTop = viewport.scrollTop;
    const viewH = viewport.clientHeight;
    const first = Math.max(0, Math.floor(scrollTop / BLOCK_H) - BUFFER);
    const last = Math.min(numBlocks - 1, Math.ceil((scrollTop + viewH) / BLOCK_H) + BUFFER);

    // remove out-of-range
    for (const key in renderedBlocks) {
        const bi = parseInt(key);
        if (bi < first || bi > last) {
            renderedBlocks[bi].remove();
            delete renderedBlocks[bi];
        }
    }

    // add missing
    const rec = SEQS[curIdx];
    const frag = document.createDocumentFragment();
    for (let bi = first; bi <= last; bi++) {
        if (renderedBlocks[bi]) continue;
        const el = createBlockEl(rec, bi);
        frag.appendChild(el);
        renderedBlocks[bi] = el;
    }
    if (frag.childNodes.length) spacer.appendChild(frag);

    // position indicator
    const midPos = Math.floor((scrollTop + viewH / 2) / BLOCK_H) * W + 1;
    const n = SEQS[curIdx].seq.length;
    posInd.textContent = midPos.toLocaleString() + ' / ' + n.toLocaleString() + ' bp';
}

// ── re-render only currently visible blocks (after search update) ──
function rerenderVisible() {
    const rec = SEQS[curIdx];
    for (const key in renderedBlocks) {
        const bi = parseInt(key);
        const newEl = createBlockEl(rec, bi);
        renderedBlocks[bi].replaceWith(newEl);
        renderedBlocks[bi] = newEl;
    }
}

// ── switch sequence ──
function selectSeq(idx) {
    curIdx = idx;
    selStart = -1; selEnd = -1;
    const rec = SEQS[curIdx];
    numBlocks = Math.ceil(rec.seq.length / W);

    buildBlastHitsCur();
    buildBlastCls(rec);
    updateBlastInfo();
    searchMarks = new Uint8Array(rec.seq.length);
    matches = [];
    activeMatch = -1;

    // reset scroll
    for (const key in renderedBlocks) { renderedBlocks[key].remove(); }
    renderedBlocks = {};
    spacer.style.height = (numBlocks * BLOCK_H) + 'px';
    viewport.scrollTop = 0;

    updateVisible();

    // re-apply search if active
    if (lastQuery) doSearch(lastQuery);
    else updateMatchInfo();
}

// ── search ──
function doSearch(query) {
    query = query.toUpperCase().replace(/[^ACGTNRYSWKMBDHV]/g, '');
    if (!query) { clearSearch(); return; }

    const rec = SEQS[curIdx];
    const seq = rec.seq;
    searchMarks = new Uint8Array(seq.length);
    matches = [];

    findAll(seq, query, matches);
    if (rcCheck.checked) {
        const rcq = rc(query);
        if (rcq !== query) findAll(seq, rcq, matches);
    }
    matches.sort((a, b) => a.pos - b.pos);

    for (const m of matches) {
        for (let i = m.pos; i < m.pos + m.len; i++) searchMarks[i] = 1;
    }

    activeMatch = matches.length > 0 ? 0 : -1;
    if (activeMatch >= 0) setActive(activeMatch);

    rerenderVisible();
    if (activeMatch >= 0) scrollToMatch(activeMatch);
    updateMatchInfo();
    lastQuery = query;
}

function findAll(seq, pat, out) {
    const plen = pat.length;
    let i = 0;
    while (true) {
        i = seq.indexOf(pat, i);
        if (i === -1) break;
        out.push({pos: i, len: plen});
        i++;
    }
}

function clearSearch() {
    const rec = SEQS[curIdx];
    searchMarks = new Uint8Array(rec.seq.length);
    matches = [];
    activeMatch = -1;
    lastQuery = '';
    rerenderVisible();
    updateMatchInfo();
}

function setActive(idx) {
    // clear previous active marks
    for (let i = 0; i < searchMarks.length; i++) {
        if (searchMarks[i] === 2) searchMarks[i] = 1;
    }
    if (idx >= 0 && idx < matches.length) {
        const m = matches[idx];
        for (let i = m.pos; i < m.pos + m.len; i++) searchMarks[i] = 2;
    }
}

function scrollToMatch(idx) {
    const m = matches[idx];
    const bi = Math.floor(m.pos / W);
    const targetTop = bi * BLOCK_H - viewport.clientHeight / 2 + BLOCK_H / 2;
    viewport.scrollTo({top: Math.max(0, targetTop), behavior: 'smooth'});
    // after scroll completes, re-render to show active highlight
    setTimeout(() => { updateVisible(); rerenderVisible(); }, 350);
}

function goNext() {
    if (!matches.length) return;
    activeMatch = (activeMatch + 1) % matches.length;
    setActive(activeMatch);
    rerenderVisible();
    scrollToMatch(activeMatch);
    updateMatchInfo();
}
function goPrev() {
    if (!matches.length) return;
    activeMatch = (activeMatch - 1 + matches.length) % matches.length;
    setActive(activeMatch);
    rerenderVisible();
    scrollToMatch(activeMatch);
    updateMatchInfo();
}

function updateMatchInfo() {
    if (!matches.length) {
        matchInfo.textContent = lastQuery ? '0 matches' : '';
        prevBtn.disabled = true;
        nextBtn.disabled = true;
    } else {
        matchInfo.textContent = (activeMatch + 1) + ' / ' + matches.length;
        prevBtn.disabled = false;
        nextBtn.disabled = false;
    }
}

// ── events ──
let searchTimer = null;
searchIn.addEventListener('input', () => {
    clearTimeout(searchTimer);
    searchTimer = setTimeout(() => doSearch(searchIn.value), 250);
});
searchIn.addEventListener('keydown', (e) => {
    if (e.key === 'Enter') { e.preventDefault(); e.shiftKey ? goPrev() : goNext(); }
    if (e.key === 'Escape') { searchIn.value = ''; clearSearch(); }
});
prevBtn.addEventListener('click', goPrev);
nextBtn.addEventListener('click', goNext);
rcCheck.addEventListener('change', () => { if (lastQuery) doSearch(lastQuery); });

seqSel.addEventListener('change', () => { selectSeq(seqSel.selectedIndex); });

gotoIn.addEventListener('keydown', (e) => {
    if (e.key === 'Enter') {
        const pos = parseInt(gotoIn.value, 10);
        if (isNaN(pos) || pos < 1) return;
        const bi = Math.floor((pos - 1) / W);
        const targetTop = bi * BLOCK_H - viewport.clientHeight / 2;
        viewport.scrollTo({top: Math.max(0, targetTop), behavior: 'smooth'});
    }
});

viewport.addEventListener('scroll', updateVisible);

document.addEventListener('keydown', (e) => {
    if ((e.ctrlKey || e.metaKey) && e.key === 'f') {
        e.preventDefault(); searchIn.focus(); searchIn.select();
    }
});

// ── BLAST from GUI ──
const blastIn     = document.getElementById('blastInput');
const blastBtn    = document.getElementById('blastBtn');
const blastPrevBtn= document.getElementById('blastPrevBtn');
const blastNextBtn= document.getElementById('blastNextBtn');
const blastClrBtn = document.getElementById('blastClearBtn');
const blastInfoEl = document.getElementById('blastInfo');

// settings panel
const blastSettingsBtn = document.getElementById('blastSettingsBtn');
const blastSettingsEl  = document.getElementById('blastSettings');
const filterMinId  = document.getElementById('filterMinId');
const filterMinPct = document.getElementById('filterMinPct');
const filterMinLen = document.getElementById('filterMinLen');
const filteredCountEl = document.getElementById('filteredCount');

let rawBlastHits = [];   // unfiltered hits from server
let blastQueryLen = 0;   // length of the query sequence (for % query filter)

blastSettingsBtn.addEventListener('click', () => {
    blastSettingsEl.classList.toggle('open');
    blastSettingsBtn.classList.toggle('active');
});

function onFilterChange() { if (rawBlastHits.length) applyFiltersAndDisplay(); }
filterMinId.addEventListener('input', onFilterChange);
filterMinPct.addEventListener('input', onFilterChange);
filterMinLen.addEventListener('input', onFilterChange);

if (blastBtn) {
    blastBtn.addEventListener('click', runBlast);
    blastIn.addEventListener('keydown', (e) => { if (e.key === 'Enter') runBlast(); });
    blastClrBtn.addEventListener('click', clearBlastHits);
    blastPrevBtn.addEventListener('click', goPrevBlast);
    blastNextBtn.addEventListener('click', goNextBlast);
}

async function runBlast() {
    let query = blastIn.value.trim();
    if (!query) return;
    // Strip FASTA header lines if pasted
    query = query.split('\n').filter(l => !l.startsWith('>')).join('');
    blastQueryLen = query.length;
    blastBtn.disabled = true;
    blastBtn.textContent = 'Running...';
    blastInfoEl.textContent = '';
    try {
        const resp = await fetch('/blast', {
            method: 'POST',
            headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({query: query})
        });
        if (!resp.ok) throw new Error('Server error ' + resp.status);
        const data = await resp.json();
        if (data.error) throw new Error(data.error);
        rawBlastHits = data.hits;
        applyFiltersAndDisplay(true);
    } catch(e) {
        blastInfoEl.textContent = 'Error: ' + e.message;
    }
    blastBtn.disabled = false;
    blastBtn.textContent = 'BLAST';
}

function filterHits(hits) {
    const minId  = parseFloat(filterMinId.value)  || 0;
    const minPct = parseFloat(filterMinPct.value)  || 0;
    const minLen = parseInt(filterMinLen.value)     || 0;
    const minLenFromPct = blastQueryLen > 0 ? blastQueryLen * minPct / 100 : 0;

    return hits.filter(h =>
        h.pi >= minId &&
        h.l >= minLen &&
        h.l >= minLenFromPct
    );
}

function applyFiltersAndDisplay(jumpToFirst) {
    const filtered = filterHits(rawBlastHits);

    filteredCountEl.textContent = filtered.length + ' / ' + rawBlastHits.length + ' hits pass filters';

    const bySeq = {};
    for (const h of filtered) {
        if (!bySeq[h.id]) bySeq[h.id] = [];
        bySeq[h.id].push({
            s: h.s, e: h.e, st: h.st, pi: h.pi, l: h.l, ev: h.ev,
            tip: 'BLAST ' + h.st + ' ' + (h.s+1) + '-' + h.e +
                 ' (' + h.pi + '%, len ' + h.l + ', e=' + h.ev.toExponential(1) + ')'
        });
    }
    BLAST_MARKS.length = 0;
    for (const id in bySeq) {
        BLAST_MARKS.push({id: id, ranges: bySeq[id]});
    }

    // If current seq has no hits, switch to first seq that does
    if (jumpToFirst) {
        const curId = SEQS[curIdx].id;
        if (!bySeq[curId] && BLAST_MARKS.length > 0) {
            const targetId = BLAST_MARKS[0].id;
            const targetIdx = SEQS.findIndex(s => s.id === targetId);
            if (targetIdx >= 0) {
                seqSel.selectedIndex = targetIdx;
                selectSeq(targetIdx);
                scrollToBlastHit(0);
                return;
            }
        }
    }

    buildBlastHitsCur();
    buildBlastCls(SEQS[curIdx]);
    rerenderVisible();
    updateBlastInfo();
    if (jumpToFirst && activeBlastHit >= 0) scrollToBlastHit(activeBlastHit);
}

function clearBlastHits() {
    BLAST_MARKS.length = 0;
    rawBlastHits = [];
    blastQueryLen = 0;
    blastHitsCur = [];
    activeBlastHit = -1;
    buildBlastCls(SEQS[curIdx]);
    rerenderVisible();
    updateBlastInfo();
    filteredCountEl.textContent = '';
    blastIn.value = '';
}

function scrollToBlastHit(idx) {
    if (idx < 0 || idx >= blastHitsCur.length) return;
    activeBlastHit = idx;
    buildBlastCls(SEQS[curIdx]);
    rerenderVisible();
    updateBlastInfo();
    const h = blastHitsCur[idx];
    const bi = Math.floor(h.s / W);
    const targetTop = bi * BLOCK_H - viewport.clientHeight / 2 + BLOCK_H / 2;
    viewport.scrollTo({top: Math.max(0, targetTop), behavior: 'smooth'});
    setTimeout(() => { updateVisible(); rerenderVisible(); }, 350);
}

function goNextBlast() {
    if (!blastHitsCur.length) return;
    activeBlastHit = (activeBlastHit + 1) % blastHitsCur.length;
    scrollToBlastHit(activeBlastHit);
}
function goPrevBlast() {
    if (!blastHitsCur.length) return;
    activeBlastHit = (activeBlastHit - 1 + blastHitsCur.length) % blastHitsCur.length;
    scrollToBlastHit(activeBlastHit);
}

function updateBlastInfo() {
    const total = BLAST_MARKS.reduce((s, b) => s + b.ranges.length, 0);
    if (!total) {
        blastInfoEl.textContent = '';
        blastPrevBtn.disabled = true;
        blastNextBtn.disabled = true;
        return;
    }
    const onSeq = blastHitsCur.length;
    if (onSeq > 0) {
        const h = blastHitsCur[activeBlastHit] || blastHitsCur[0];
        blastInfoEl.textContent = (activeBlastHit + 1) + ' / ' + onSeq +
            ' on ' + SEQS[curIdx].id + ' (' + total + ' total)' +
            '  |  ' + h.pi + '% id, len ' + h.l + ', ' + h.st +
            ' strand, pos ' + (h.s+1) + '-' + h.e + ', e=' + h.ev.toExponential(1);
        blastPrevBtn.disabled = false;
        blastNextBtn.disabled = false;
    } else {
        blastInfoEl.textContent = '0 on ' + SEQS[curIdx].id + ' (' + total + ' total)';
        blastPrevBtn.disabled = true;
        blastNextBtn.disabled = true;
    }
}

// ── sequence selection (SnapGene-style) ──
let selStart = -1;
let selEnd = -1;
let selecting = false;
const posTip = document.getElementById('posTip');

function getSeqPos(e) {
    const el = e.target.closest('[data-p]');
    if (el) return parseInt(el.dataset.p);
    // Fallback: estimate from Y position
    const y = e.clientY - spacer.getBoundingClientRect().top + viewport.scrollTop;
    const bi = Math.floor(y / BLOCK_H);
    const n = SEQS[curIdx].seq.length;
    return Math.max(0, Math.min(bi * W + W - 1, n - 1));
}

function showTip(e, text) {
    posTip.textContent = text;
    posTip.style.display = 'block';
    posTip.style.left = (e.clientX + 14) + 'px';
    posTip.style.top = (e.clientY - 28) + 'px';
}
function hideTip() { posTip.style.display = 'none'; }

spacer.addEventListener('mousedown', (e) => {
    if (e.button !== 0) return;
    const pos = getSeqPos(e);
    if (pos < 0) return;
    selStart = pos;
    selEnd = pos;
    selecting = true;
    rerenderVisible();
    showTip(e, '' + (pos + 1));
    e.preventDefault();
});

document.addEventListener('mousemove', (e) => {
    if (selecting) {
        const pos = getSeqPos(e);
        if (pos < 0) return;
        if (pos !== selEnd) {
            selEnd = pos;
            rerenderVisible();
        }
        const lo = Math.min(selStart, selEnd) + 1;
        const hi = Math.max(selStart, selEnd) + 1;
        const len = hi - lo + 1;
        if (lo === hi) showTip(e, '' + lo);
        else showTip(e, lo + ' .. ' + hi + '  =  ' + len + ' bp');
        return;
    }
    // Hover: show position of base under cursor
    const el = e.target.closest('[data-p]');
    if (el) {
        showTip(e, '' + (parseInt(el.dataset.p) + 1));
    } else {
        hideTip();
    }
});

document.addEventListener('mouseup', () => {
    if (selecting) {
        selecting = false;
        // Keep tip visible showing final selection range
    }
});

// Hide tip when leaving the sequence area
viewport.addEventListener('mouseleave', hideTip);

// Clear selection when clicking toolbar inputs
document.getElementById('toolbar').addEventListener('mousedown', () => {
    if (selStart >= 0) { selStart = -1; selEnd = -1; rerenderVisible(); }
    hideTip();
});

// Copy forward strand only
document.addEventListener('copy', (e) => {
    if (selStart < 0 || selStart === selEnd) return;
    e.preventDefault();
    const lo = Math.min(selStart, selEnd);
    const hi = Math.max(selStart, selEnd);
    const text = SEQS[curIdx].seq.substring(lo, hi + 1);
    e.clipboardData.setData('text/plain', text);
});

// ── init ──
measureBlockH();
selectSeq(0);
"""


def make_server_handler(html_content, fasta_path):
    """Create an HTTP request handler with access to the HTML and FASTA path."""
    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            if self.path == "/" or self.path == "":
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.end_headers()
                self.wfile.write(html_content.encode("utf-8"))
            else:
                self.send_error(404)

        def do_POST(self):
            if self.path == "/blast":
                try:
                    length = int(self.headers.get("Content-Length", 0))
                    body = json.loads(self.rfile.read(length))
                    query = body.get("query", "").strip()
                    if not query:
                        raise ValueError("Empty query")
                    evalue = float(body.get("evalue", 10.0))
                    word_size = int(body.get("wordSize", 7))
                    hits = run_blast(str(fasta_path), query,
                                     evalue=evalue, word_size=word_size)
                    resp = json.dumps({"hits": hits})
                    self.send_response(200)
                except Exception as e:
                    resp = json.dumps({"error": str(e)})
                    self.send_response(500)
                self.send_header("Content-Type", "application/json")
                self.end_headers()
                self.wfile.write(resp.encode())
            else:
                self.send_error(404)

        def handle_one_request(self):
            try:
                super().handle_one_request()
            except (ConnectionResetError, BrokenPipeError):
                pass

        def log_message(self, fmt, *a):
            msg = str(a[0]) if a else ""
            if "/blast" in msg:
                super().log_message(fmt, *a)

    return Handler


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", help="Reference FASTA file to visualize")
    ap.add_argument("-o", "--output", default=None,
                    help="Output HTML file (default: <fasta stem>.html)")
    ap.add_argument("-w", "--width", type=int, default=100,
                    help="Bases per line (default: 100)")
    ap.add_argument("--serve", action="store_true",
                    help="Start a local server with live BLAST support.")
    ap.add_argument("--port", type=int, default=8050,
                    help="Port for --serve mode (default: 8050)")
    ap.add_argument("--blast", default=None,
                    help="BLAST query: FASTA path or literal sequence. Hits shown green.")
    ap.add_argument("--evalue", type=float, default=10.0)
    ap.add_argument("--word-size", type=int, default=7)
    args = ap.parse_args()

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        sys.exit(f"File not found: {fasta_path}")

    records = parse_fasta(fasta_path)

    # BLAST
    blast_marks = {}
    blast_hits_all = []
    if args.blast:
        blast_hits_all = run_blast(str(fasta_path), args.blast,
                                    evalue=args.evalue, word_size=args.word_size)
        for h in blast_hits_all:
            blast_marks.setdefault(h["id"], []).append({
                "s": h["s"], "e": h["e"], "st": h["st"],
                "tip": f"BLAST {h['st']} {h['s']+1}-{h['e']} "
                       f"({h['pi']}%, len {h['l']}, e={h['ev']:.1e})",
            })

    # Build JS data
    seq_data = []
    blast_data = []
    for name, seq in records:
        seq_data.append({
            "id": name, "name": name,
            "seq": seq, "comp": complement_str(seq),
        })
        bm = blast_marks.get(name, [])
        if bm:
            blast_data.append({"id": name, "ranges": bm})

    # Sequence selector options
    options_html = "".join(
        f'<option value="{i}">{html.escape(r["name"])} ({len(r["seq"]):,} bp)</option>'
        for i, r in enumerate(seq_data)
    )

    # Blast summary
    blast_summary = ""
    if args.blast:
        q_disp = args.blast if len(args.blast) < 60 else args.blast[:57] + "..."
        if blast_hits_all:
            rows = "".join(
                f"<tr><td>{html.escape(h['id'])}</td><td>{h['st']}</td>"
                f"<td>{h['s']+1}-{h['e']}</td><td>{h['pi']}%</td>"
                f"<td>{h['l']}</td><td>{h['ev']:.1e}</td></tr>"
                for h in blast_hits_all
            )
            table = (
                "<table><tr><th>subject</th><th>strand</th><th>range</th>"
                f"<th>pident</th><th>len</th><th>evalue</th></tr>{rows}</table>"
            )
        else:
            table = "<em>no hits</em>"
        blast_summary = (
            f'<div class="blast-summary">BLAST query: '
            f'<code>{html.escape(q_disp)}</code> &mdash; '
            f'{len(blast_hits_all)} hit(s)<br>{table}</div>'
        )

    # Legend
    legend_parts = ['<span class="swatch hl-search">AAAA</span> search match']
    if args.blast:
        legend_parts.append('<span class="swatch hl-blast">AAAA</span> BLAST hit')
    legend_parts.append('<span class="swatch hl-active">AAAA</span> active match')
    legend_html = '<div class="legend">' + ' &nbsp; '.join(legend_parts) + '</div>'

    html_doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{html.escape(fasta_path.name)} — Sequence Viewer</title>
<style>{CSS}</style>
</head>
<body>

<div id="toolbar">
  <label for="seqSelect">Sequence</label>
  <select id="seqSelect">{options_html}</select>

  <div class="sep"></div>

  <label for="searchInput">Search</label>
  <input type="text" id="searchInput" placeholder="Enter DNA sequence..." spellcheck="false" autocomplete="off">
  <label id="rcLabel"><input type="checkbox" id="rcCheck" checked> + RC</label>
  <button id="prevBtn" disabled title="Previous match (Shift+Enter)">&laquo; Prev</button>
  <button id="nextBtn" disabled title="Next match (Enter)">Next &raquo;</button>
  <span id="matchInfo"></span>

  <div class="sep"></div>

  <label for="gotoInput">Go to</label>
  <input type="text" id="gotoInput" placeholder="position" spellcheck="false">

  <div class="sep"></div>

  <label for="blastInput">BLAST</label>
  <input type="text" id="blastInput" placeholder="Paste sequence to BLAST..." spellcheck="false" autocomplete="off">
  <button id="blastBtn">BLAST</button>
  <button id="blastPrevBtn" disabled title="Previous BLAST hit">&laquo;</button>
  <button id="blastNextBtn" disabled title="Next BLAST hit">&raquo;</button>
  <button id="blastClearBtn" title="Clear BLAST highlights">Clear</button>
  <button id="blastSettingsBtn" title="BLAST filter settings">&#9881;</button>
  <span id="blastInfo"></span>

  <span id="posIndicator"></span>
</div>

<div id="blastSettings">
  <span class="filter-title">BLAST Filters:</span>
  <label>Min % identity <input type="number" id="filterMinId" value="0" min="0" max="100" step="1"></label>
  <label>Min % query length <input type="number" id="filterMinPct" value="0" min="0" max="100" step="5"></label>
  <label>Min alignment bp <input type="number" id="filterMinLen" value="0" min="0" step="10"></label>
  <span id="filteredCount"></span>
</div>

<div id="infobar">
{legend_html}
{blast_summary}
</div>

<div id="viewport">
  <div id="spacer"></div>
</div>
<div id="posTip"></div>

<script>
const SEQS = {json.dumps(seq_data, separators=(',', ':'))};
const BLAST_MARKS = {json.dumps(blast_data, separators=(',', ':'))};
const LINE_WIDTH = {args.width};
{JS}
</script>
</body>
</html>
"""

    if args.serve:
        handler = make_server_handler(html_doc, fasta_path)
        server = HTTPServer(("127.0.0.1", args.port), handler)
        url = f"http://127.0.0.1:{args.port}"
        print(f"Serving {fasta_path.name} at {url}")
        print(f"  {len(records)} sequence(s), "
              f"{sum(len(s) for _, s in records):,} bp total")
        print(f"  BLAST endpoint: POST {url}/blast")
        print("  Press Ctrl+C to stop.")
        import subprocess as _sp
        try:
            _sp.Popen(
                ["python3", "-c", f"import webbrowser; webbrowser.open('{url}')"],
                stdout=_sp.DEVNULL, stderr=_sp.DEVNULL,
            )
        except Exception:
            pass
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            print("\nStopped.")
            server.server_close()
    else:
        out_path = Path(args.output) if args.output else fasta_path.with_suffix(".html")
        out_path.write_text(html_doc)
        print(f"Wrote {out_path} ({len(records)} sequence(s), "
              f"{sum(len(s) for _, s in records):,} bp total)")
        if args.blast:
            print(f"  BLAST hits: {len(blast_hits_all)}")


if __name__ == "__main__":
    main()
