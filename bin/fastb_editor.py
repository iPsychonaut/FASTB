#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FASTB Editor (PyQt)

Minimal editor that treats a .fastb file like text:
- Open a .fastb -> decode to FASTA-like text in the editor.
- Save -> re-encode to .fastb with a .bak backup and integrity validation.

Dependencies:
    - PyQt5
    - bitarray
    - Your existing fastb.py module in the same folder (providing:
        fastb_transform, read_fastb_file, save_bitarray_to_file)

Usage:
    python fastb_editor.py                # launches empty, use File > Open
    python fastb_editor.py /path/file.fastb

Editable text format (FASTB-EDIT-TXT v1):
    - Comments start with '#'
    - Header lines start with '>' and may include an advisory NUC hint:
        >DESCRIPTION TEXT | NUC=DNA
        >Another record without hint
    - Sequence lines follow (A/C/G/T/U/N only), any wrapping ok
    - Blank line separates records

Round-trip is validated on Save by decoding the newly written .fastb and
ensuring at least one record exists (and basic sanity of lengths).
"""

import sys
import os
import traceback
from pathlib import Path
from typing import List, Tuple, Optional

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QFileDialog, QMessageBox,
    QAction, QPlainTextEdit, QWidget, QVBoxLayout, QLabel, QStatusBar
)
from PyQt5.QtCore import Qt

# Import your codec module
# Ensure fastb.py is in the same directory or on PYTHONPATH
import importlib
fastb = importlib.import_module("fastb")  # your existing module

WRAP = 80  # purely cosmetic when exporting text; editor itself doesn't hard-wrap


def decode_fastb_to_text(fastb_path: Path) -> str:
    """
    Decode a .fastb file into FASTB-EDIT-TXT v1 text.

    Args:
        fastb_path: Path to .fastb file.

    Returns:
        str: Text content for the editor.
    """
    records = fastb.read_fastb_file(str(fastb_path))
    lines = []
    lines.append("# FASTB-EDIT-TXT v1")
    lines.append("# Lines starting with '>' are headers. Optional meta: '| NUC=DNA|RNA'")
    lines.append("# Blank line separates records. Sequences are raw A/C/G/T/U/N.\n")
    for desc, seq, nuc in records:
        header = f">{desc}"
        if nuc in {"DNA", "RNA"}:
            header += f" | NUC={nuc}"
        lines.append(header)
        for i in range(0, len(seq), WRAP):
            lines.append(seq[i:i + WRAP])
        lines.append("")  # blank line between records
    return "\n".join(lines).rstrip() + "\n"


def _parse_edit_text(text: str) -> List[Tuple[str, str]]:
    """
    Parse FASTB-EDIT-TXT v1 text into a list of (description, sequence).

    Notes:
        - NUC hint is advisory; encoding deduces T/U itself.
        - Comments (#) ignored. Blank line flushes the current record.

    Returns:
        list of (desc, sequence)
    """
    records: List[Tuple[str, str]] = []
    desc: Optional[str] = None
    seq_chunks: List[str] = []

    def flush():
        nonlocal desc, seq_chunks
        if desc is None:
            return
        seq = "".join(seq_chunks).replace(" ", "").replace("\t", "")
        if not seq:
            raise ValueError(f"Empty sequence for record: {desc}")
        records.append((desc, seq))
        desc = None
        seq_chunks = []

    for raw in text.splitlines():
        line = raw.rstrip("\n")
        if not line:
            flush()
            continue
        if line.startswith("#"):
            continue
        if line.startswith(">"):
            flush()
            header = line[1:].strip()
            if " | " in header:
                parts = header.split(" | ", 1)
                desc = parts[0].strip()
                # meta is ignored for encoding; left for humans
            else:
                desc = header
        else:
            seq_chunks.append(line.strip())

    flush()
    if not records:
        raise ValueError("No records parsed from editor text.")
    return records


def encode_text_to_fastb_bytes(text: str) -> bytes:
    """
    Convert FASTB-EDIT-TXT text into aggregated FASTB bytes.

    Returns:
        bytes: Aggregated FASTB data suitable for writing to .fastb
    """
    pairs = _parse_edit_text(text)  # [(desc, seq), ...]
    # Build aggregate bitarray using your codec
    from bitarray import bitarray
    agg = bitarray()
    for desc, seq in pairs:
        rec_bits = fastb.fastb_transform(seq, "encode", desc)
        agg.extend(rec_bits)
    return agg.tobytes()


def validate_fastb_roundtrip(fastb_path: Path) -> None:
    """
    Decode newly written .fastb to ensure it's structurally valid.

    Raises:
        ValueError if file fails basic sanity checks.
    """
    recs = fastb.read_fastb_file(str(fastb_path))
    if not recs:
        raise ValueError("Validation failed: no records decoded.")
    # trivial sanity: descriptions non-empty, lengths >= 1
    for i, (desc, seq, nuc) in enumerate(recs, start=1):
        if not isinstance(desc, str) or desc == "":
            raise ValueError(f"Validation failed: empty description at record {i}.")
        if not isinstance(seq, str) or len(seq) == 0:
            raise ValueError(f"Validation failed: empty sequence at record {i}.")
        if nuc not in {"DNA", "RNA"}:
            raise ValueError(f"Validation failed: invalid nucleotide type at record {i}: {nuc}")


class FastbEditor(QMainWindow):
    """
    Tiny PyQt editor that decodes/encodes .fastb transparently.
    """

    def __init__(self, path: Optional[Path] = None):
        super().__init__()
        self.setWindowTitle("FASTB Editor")
        self.resize(900, 700)

        self._fastb_path: Optional[Path] = None  # original .fastb path

        # Widgets
        self.editor = QPlainTextEdit()
        self.editor.setLineWrapMode(QPlainTextEdit.NoWrap)
        self.editor.setTabChangesFocus(False)
        self.editor.textChanged.connect(self._on_text_changed)

        self.status = QStatusBar()
        self.setStatusBar(self.status)
        self._dirty = False

        container = QWidget()
        layout = QVBoxLayout(container)
        layout.addWidget(self.editor)
        self.setCentralWidget(container)

        self._build_menu()

        if path:
            self.open_fastb(path)

    # ---- Menu ----
    def _build_menu(self):
        m_file = self.menuBar().addMenu("&File")

        act_open = QAction("&Open .fastb…", self)
        act_open.triggered.connect(self.action_open)
        m_file.addAction(act_open)

        act_save = QAction("&Save", self)
        act_save.setShortcut("Ctrl+S")
        act_save.triggered.connect(self.action_save)
        m_file.addAction(act_save)

        act_save_as = QAction("Save &As…", self)
        act_save_as.triggered.connect(self.action_save_as)
        m_file.addAction(act_save_as)

        m_file.addSeparator()

        act_export = QAction("Export text…", self)
        act_export.triggered.connect(self.action_export_text)
        m_file.addAction(act_export)

        m_file.addSeparator()

        act_quit = QAction("&Quit", self)
        act_quit.setShortcut("Ctrl+Q")
        act_quit.triggered.connect(self.close)
        m_file.addAction(act_quit)

        m_tools = self.menuBar().addMenu("&Tools")

        act_validate = QAction("&Validate current text encodes", self)
        act_validate.triggered.connect(self.action_validate_text)
        m_tools.addAction(act_validate)

        act_reload = QAction("&Reload from .fastb", self)
        act_reload.triggered.connect(self.action_reload_from_disk)
        m_tools.addAction(act_reload)

        m_view = self.menuBar().addMenu("&View")

        act_wrap = QAction("Toggle soft wrap", self, checkable=True)
        act_wrap.triggered.connect(self.action_toggle_wrap)
        m_view.addAction(act_wrap)

    # ---- Actions ----
    def action_open(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open .fastb", "", "FASTB (*.fastb)")
        if not path:
            return
        self.open_fastb(Path(path))

    def action_save(self):
        if self._fastb_path is None:
            self.action_save_as()
        else:
            self.save_to_fastb(self._fastb_path)

    def action_save_as(self):
        if self._fastb_path is None:
            default = os.path.expanduser("~/untitled.fastb")
        else:
            default = str(self._fastb_path)
        path, _ = QFileDialog.getSaveFileName(self, "Save .fastb", default, "FASTB (*.fastb)")
        if not path:
            return
        self.save_to_fastb(Path(path))

    def action_export_text(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export decoded text", "", "Text (*.txt)")
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(self.editor.toPlainText())
            self.status.showMessage(f"Exported text: {path}", 4000)
        except Exception as e:
            self._error("Failed to export text", e)

    def action_validate_text(self):
        try:
            _ = encode_text_to_fastb_bytes(self.editor.toPlainText())
            QMessageBox.information(self, "Validate", "Text encodes successfully.")
        except Exception as e:
            self._error("Validation failed while encoding text", e)

    def action_reload_from_disk(self):
        if self._fastb_path is None:
            QMessageBox.information(self, "Reload", "No .fastb file is associated yet.")
            return
        self.open_fastb(self._fastb_path)

    def action_toggle_wrap(self, checked: bool):
        self.editor.setLineWrapMode(QPlainTextEdit.WidgetWidth if checked else QPlainTextEdit.NoWrap)

    # ---- Core IO ----
    def open_fastb(self, path: Path):
        try:
            text = decode_fastb_to_text(path)
        except Exception as e:
            self._error(f"Failed to decode: {path}", e)
            return
        self._fastb_path = path
        self._set_text_clean(text)
        self._update_title()
        self.status.showMessage(f"Opened {path}", 4000)

    def save_to_fastb(self, path: Path):
        # Encode editor text to bytes, write with backup, validate
        try:
            data = encode_text_to_fastb_bytes(self.editor.toPlainText())
        except Exception as e:
            self._error("Encoding failed. Save aborted.", e)
            return

        try:
            # backup if exists
            if path.exists():
                backup = path.with_suffix(path.suffix + ".bak")
                backup.write_bytes(path.read_bytes())
            # write via your helper for consistency
            from bitarray import bitarray
            ba = bitarray()
            ba.frombytes(data)
            fastb.save_bitarray_to_file(ba, str(path))
        except Exception as e:
            self._error(f"Failed to write .fastb: {path}", e)
            return

        # Validate by decoding what we wrote
        try:
            validate_fastb_roundtrip(path)
        except Exception as e:
            self._error(f"Wrote file but validation failed: {path}", e)
            return

        self._fastb_path = path
        self._dirty = False
        self._update_title()
        self.status.showMessage(f"Saved {path}", 4000)

    # ---- Helpers ----
    def _on_text_changed(self):
        self._dirty = True
        self._update_title()

    def _set_text_clean(self, text: str):
        self.editor.blockSignals(True)
        self.editor.setPlainText(text)
        self.editor.blockSignals(False)
        self._dirty = False

    def _update_title(self):
        name = str(self._fastb_path) if self._fastb_path else "(unsaved)"
        mark = "*" if self._dirty else ""
        self.setWindowTitle(f"FASTB Editor - {name}{mark}")

    def _error(self, title: str, exc: Exception):
        tb = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
        QMessageBox.critical(self, title, f"{exc}\n\nDetails:\n{tb}")


def main():
    app = QApplication(sys.argv)
    path = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else None
    w = FastbEditor(path if path and path.suffix == ".fastb" else None)
    w.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
