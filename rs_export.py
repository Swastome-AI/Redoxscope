"""
RedoxScope - plot styling, publication-quality figure export and PDF report.
No Streamlit code in this file.
"""
import datetime
import html
import importlib.util
import io
import os
import re
import textwrap

import plotly.graph_objects as go
import plotly.io as pio

# ----------------------------------------------------------------------------
# Branding / footer
# ----------------------------------------------------------------------------
FOOTER_1 = ("RedoxScope is actively maintained and developed by Swastome-AI "
            "under Swastome Private Limited.")
FOOTER_2 = ("© 2026 Swastome Private Limited. All rights reserved. Swastome and the Swastome logo are "
            "the property of Swastome Private Limited. No part of this app may be copied, modified or "
            "distributed without written permission. Email to management@swastome.com for "
            "collaboration and licensing enquiries.")

# ----------------------------------------------------------------------------
# Plot style: light mode only, black text, big fonts, grid lines
# ----------------------------------------------------------------------------
BLACK = "#000000"
GRID = "#BFBFBF"
FONT = "Arial, Helvetica, sans-serif"
# Colour-blind-safe palette (Okabe-Ito, without yellow)
COLORWAY = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#999999", "#000000"]

pio.templates["redoxscope"] = go.layout.Template(layout=dict(colorway=COLORWAY))
pio.templates.default = "plotly_white+redoxscope"

# Light colour scales: dark text stays readable on every cell
SCALE_REDS = ["#FFFFFF", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A"]
SCALE_ORANGES = ["#FFFFFF", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C"]
SCALE_BLUES_REV = ["#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7", "#F7FBFF"]

PUB_W, PUB_H = 1000, 700  # export size in pixels; PNG/TIFF at scale 3 = 3000 px wide = 300 dpi at 10 in
FORMATS = {
    "PNG – 300 dpi (slides, web)": ("png", "image/png"),
    "TIFF – 300 dpi (journals)": ("tiff", "image/tiff"),
    "SVG – vector (editable)": ("svg", "image/svg+xml"),
    "PDF – vector (print)": ("pdf", "application/pdf"),
}


def _num(values):
    out = []
    for v in values:
        try:
            f = float(v)
        except (TypeError, ValueError):
            continue
        if f == f:
            out.append(f)
    return out


def _bar_headroom(fig):
    """Give bars with value labels room above/below so the labels are never cut off."""
    bars = [t for t in fig.data if t.type == "bar" and t.text is not None and t.y is not None]
    if not bars:
        return
    ys = _num(v for t in bars for v in t.y)
    if not ys:
        return
    lo, hi = min(0.0, min(ys)), max(0.0, max(ys))
    pad = 0.15 * max(hi - lo, 1e-9)
    ya = fig.layout.yaxis
    if ya.range is not None:
        r = list(ya.range)
        if r[1] > r[0]:
            r[1] = max(r[1], hi + pad)
            fig.update_yaxes(range=r)
        return
    if ya.autorange == "reversed":
        fig.update_yaxes(range=[hi, lo - pad if lo < 0 else lo])
    else:
        fig.update_yaxes(range=[lo - pad if lo < 0 else 0, hi + pad if hi > 0 else 0])


def style_fig(fig):
    """Apply the RedoxScope look: white background, black text, large fonts, grid lines."""
    title = fig.layout.title.text or ""
    lines = textwrap.wrap(title, 70) if title else []
    n = max(1, len(lines))
    fig.update_layout(
        template="plotly_white+redoxscope",
        title=dict(text="<br>".join(lines), x=0.0, xanchor="left",
                   font=dict(size=21, color=BLACK, family=FONT)),
        font=dict(family=FONT, size=16, color=BLACK),
        paper_bgcolor="white", plot_bgcolor="white",
        margin=dict(t=40 + 34 * n, b=75, l=90, r=40),
        legend=dict(font=dict(size=15, color=BLACK), title=dict(font=dict(size=16, color=BLACK)),
                    bgcolor="white", bordercolor=BLACK, borderwidth=1),
        coloraxis_colorbar=dict(tickfont=dict(size=14, color=BLACK),
                                title=dict(font=dict(size=15, color=BLACK)),
                                outlinecolor=BLACK, outlinewidth=1),
    )
    axis = dict(showgrid=True, gridcolor=GRID, gridwidth=1, showline=True, linecolor=BLACK,
                linewidth=1.5, mirror=True, ticks="outside", tickcolor=BLACK, ticklen=6,
                tickfont=dict(size=15, color=BLACK), title_font=dict(size=17, color=BLACK),
                automargin=True)
    fig.update_xaxes(zeroline=False, **axis)
    fig.update_yaxes(zeroline=True, zerolinecolor="#404040", zerolinewidth=1.5, **axis)

    # heatmaps: black numbers, thin cell borders
    fig.update_traces(textfont=dict(size=15, color=BLACK), xgap=1, ygap=1,
                      selector=dict(type="heatmap"))
    # bars: black value labels outside the bar
    has_bar_text = False
    for t in fig.data:
        if t.type == "bar" and t.text is not None:
            t.textposition = "outside"
            t.textfont = dict(size=15, color=BLACK)
            t.cliponaxis = False
            has_bar_text = True
    if has_bar_text:
        if fig.layout.barmode in (None, "relative"):
            fig.update_layout(barmode="overlay")  # one bar per category, so labels sit on top
        _bar_headroom(fig)
    # scatter labels
    fig.update_traces(textfont=dict(size=15, color=BLACK),
                      selector=lambda t: t.type == "scatter" and t.mode is not None and "text" in t.mode)
    # radar
    if any(t.type == "scatterpolar" for t in fig.data):
        fig.update_polars(
            bgcolor="white",
            radialaxis=dict(gridcolor=GRID, linecolor=BLACK, tickfont=dict(size=14, color=BLACK)),
            angularaxis=dict(gridcolor=GRID, linecolor=BLACK, tickfont=dict(size=15, color=BLACK)))
    return fig


def plain_title(fig):
    return (fig.layout.title.text or "").replace("<br>", " ")


# ----------------------------------------------------------------------------
# Figure export (needs the 'kaleido' package)
# ----------------------------------------------------------------------------
def _pil_dpi(png_bytes, fmt):
    from PIL import Image
    im = Image.open(io.BytesIO(png_bytes))
    out = io.BytesIO()
    if fmt == "tiff":
        im.convert("RGB").save(out, format="TIFF", dpi=(300, 300), compression="tiff_lzw")
    else:
        im.save(out, format="PNG", dpi=(300, 300))
    return out.getvalue()


def fig_bytes(fig, fmt, scale=3):
    """Return the figure as png / tiff / svg / pdf bytes at a fixed publication size."""
    f = go.Figure(fig)
    w, h = PUB_W, int(fig.layout.height or PUB_H)
    f.update_layout(width=w, height=h)
    if fmt in ("png", "tiff"):
        png = f.to_image(format="png", width=w, height=h, scale=scale)
        try:
            return _pil_dpi(png, fmt)
        except Exception:  # noqa
            if fmt == "tiff":
                raise
            return png
    return f.to_image(format=fmt, width=w, height=h)


# ----------------------------------------------------------------------------
# PDF report (reportlab)
# ----------------------------------------------------------------------------
def safe_name(name):
    return re.sub(r"[^A-Za-z0-9]+", "_", name or "").strip("_") or "user"


def _register_fonts():
    """Use DejaVu Sans (full unicode) when available, else Helvetica."""
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.lib.fonts import addMapping
    dirs = ["/usr/share/fonts/truetype/dejavu", "/usr/share/fonts/dejavu",
            "/Library/Fonts", "C:/Windows/Fonts"]
    spec = importlib.util.find_spec("matplotlib")
    if spec and spec.origin:
        dirs.append(os.path.join(os.path.dirname(spec.origin), "mpl-data", "fonts", "ttf"))
    for d in dirs:
        reg, bold = os.path.join(d, "DejaVuSans.ttf"), os.path.join(d, "DejaVuSans-Bold.ttf")
        if os.path.exists(reg) and os.path.exists(bold):
            try:
                pdfmetrics.registerFont(TTFont("RS", reg))
                pdfmetrics.registerFont(TTFont("RS-Bold", bold))
                addMapping("RS", 0, 0, "RS")
                addMapping("RS", 1, 0, "RS-Bold")
                addMapping("RS", 0, 1, "RS")
                addMapping("RS", 1, 1, "RS-Bold")
                return "RS", "RS-Bold", True
            except Exception:  # noqa
                pass
    return "Helvetica", "Helvetica-Bold", False


_ASCII_FALLBACK = {"→": "->", "−": "-", "≤": "<=", "≥": ">=", "≈": "~", "×": "x", "…": "...",
                   "“": '"', "”": '"', "’": "'", "‘": "'", "•": "-", "✓": "Yes", "✗": "No"}
_SYMBOLS = {"✅": "Yes", "❌": "No", "🟡": "?", "⚠️": "!", "\u2009": " "}


def _clean(text, unicode_ok):
    s = str(text)
    for k, v in _SYMBOLS.items():
        s = s.replace(k, v)
    if not unicode_ok:
        for k, v in _ASCII_FALLBACK.items():
            s = s.replace(k, v)
        s = s.encode("cp1252", "replace").decode("cp1252")
    return re.sub(r"[\U00010000-\U0010FFFF]", "", s)  # drop remaining emoji


def _fmt(v):
    try:
        import pandas as pd
        if pd.isna(v):
            return ""
    except Exception:  # noqa
        pass
    if isinstance(v, float):
        return f"{v:.3f}".rstrip("0").rstrip(".")
    return str(v)


def build_pdf(user_name, app_name, app_tag, findings, overview, tables, figures, notes=None):
    """
    findings : list of strings (may contain <b>..</b>)
    overview : list of (label, value) pairs
    tables   : dict name -> DataFrame
    figures  : list of (caption, png_bytes)
    notes    : optional list of extra strings (e.g. 'plots could not be rendered')
    returns  : PDF bytes
    """
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_CENTER
    from reportlab.lib.pagesizes import A4, landscape
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.units import cm
    from reportlab.lib.utils import ImageReader
    from reportlab.platypus import (BaseDocTemplate, Frame, Image, KeepTogether, PageBreak,
                                    PageTemplate, Paragraph, Spacer, Table, TableStyle)

    reg, bold, uni = _register_fonts()
    C = lambda t: _clean(t, uni)  # noqa
    E = lambda t: html.escape(C(t))  # noqa
    # Findings use <b> tags only: escape everything, then restore them
    F = lambda t: html.escape(C(t)).replace("&lt;b&gt;", "<b>").replace("&lt;/b&gt;", "</b>")  # noqa

    navy = colors.HexColor("#1F3A5F")
    st_title = ParagraphStyle("t", fontName=bold, fontSize=26, leading=32, textColor=navy, spaceAfter=6)
    st_sub = ParagraphStyle("s", fontName=reg, fontSize=13, leading=17, textColor=colors.black, spaceAfter=4)
    st_h1 = ParagraphStyle("h1", fontName=bold, fontSize=17, leading=22, textColor=navy, spaceBefore=10, spaceAfter=8)
    st_h2 = ParagraphStyle("h2", fontName=bold, fontSize=12.5, leading=16, textColor=colors.black, spaceBefore=8, spaceAfter=4)
    st_body = ParagraphStyle("b", fontName=reg, fontSize=10.5, leading=15, textColor=colors.black, spaceAfter=4)
    st_cap = ParagraphStyle("c", fontName=reg, fontSize=10, leading=13, textColor=colors.black, alignment=TA_CENTER, spaceBefore=4)
    st_cell = ParagraphStyle("cell", fontName=reg, fontSize=8, leading=10, textColor=colors.black)
    st_head = ParagraphStyle("hd", fontName=bold, fontSize=8, leading=10, textColor=colors.white)
    st_foot = ParagraphStyle("f", fontName=reg, fontSize=6.5, leading=8, textColor=colors.black, alignment=TA_CENTER)

    page = landscape(A4)
    lm = rm = 1.6 * cm
    tm, bm = 1.5 * cm, 2.6 * cm
    avail_w = page[0] - lm - rm
    avail_h = page[1] - tm - bm
    date_s = datetime.datetime.now().strftime("%d %B %Y, %H:%M")
    who = C(user_name)

    def on_page(canv, doc):
        canv.saveState()
        canv.setStrokeColor(colors.HexColor("#999999"))
        canv.line(lm, bm - 0.25 * cm, page[0] - rm, bm - 0.25 * cm)
        canv.setFont(reg, 8)
        canv.setFillColor(colors.black)
        canv.drawString(lm, bm - 0.7 * cm, f"{app_name} report  |  Prepared by {who}  |  {date_s}")
        canv.drawRightString(page[0] - rm, bm - 0.7 * cm, f"Page {doc.page}")
        p = Paragraph(E(FOOTER_1) + "<br/>" + E(FOOTER_2), st_foot)
        w, h = p.wrap(avail_w, 1.2 * cm)
        p.drawOn(canv, lm, 0.35 * cm)
        canv.restoreState()

    buf = io.BytesIO()
    doc = BaseDocTemplate(buf, pagesize=page, leftMargin=lm, rightMargin=rm, topMargin=tm, bottomMargin=bm,
                          title=f"{app_name} report - {who}", author=who, subject=app_tag)
    doc.addPageTemplates([PageTemplate(id="p", frames=[Frame(lm, bm, avail_w, avail_h, id="f")], onPage=on_page)])

    S = []
    # ---- cover / summary
    S += [Spacer(1, 0.6 * cm), Paragraph(E(f"{app_name} Analysis Report"), st_title),
          Paragraph(E(app_tag), st_sub), Spacer(1, 0.3 * cm),
          Paragraph(f"<b>Prepared by:</b> {E(user_name)}", st_sub),
          Paragraph(f"<b>Date:</b> {E(date_s)}", st_sub), Spacer(1, 0.4 * cm)]
    if overview:
        rows = [[Paragraph(f"<b>{E(k)}</b>", st_cell), Paragraph(E(v), st_cell)] for k, v in overview]
        t = Table(rows, colWidths=[6 * cm, 10 * cm], hAlign="LEFT")
        t.setStyle(TableStyle([("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#999999")),
                               ("BACKGROUND", (0, 0), (0, -1), colors.HexColor("#EEF2F7")),
                               ("VALIGN", (0, 0), (-1, -1), "MIDDLE")]))
        S += [t, Spacer(1, 0.4 * cm)]
    S.append(Paragraph("Key findings", st_h1))
    for f in findings:
        S.append(Paragraph("• " + F(f), st_body))
    S.append(Paragraph(
        "<i>All results are computer predictions (ProTox-3 and CB-Dock2) for learning and hypothesis building. "
        "They are not laboratory evidence.</i>", st_body))
    for n in (notes or []):
        S.append(Paragraph("<b>Note:</b> " + E(n), st_body))

    # ---- figures
    if figures:
        S += [PageBreak(), Paragraph("Figures", st_h1)]
        for i, (cap, png) in enumerate(figures, 1):
            iw, ih = ImageReader(io.BytesIO(png)).getSize()
            max_w, max_h = avail_w * 0.92, avail_h - 2.4 * cm
            sc = min(max_w / iw, max_h / ih)
            block = KeepTogether([Image(io.BytesIO(png), width=iw * sc, height=ih * sc, hAlign="CENTER"),
                                  Paragraph(f"<b>Figure {i}.</b> {E(cap)}", st_cap)])
            S += [block, PageBreak() if i < len(figures) else Spacer(1, 0.2 * cm)]

    # ---- tables
    if tables:
        S += [PageBreak(), Paragraph("Tables", st_h1)]
        for ti, (name, df) in enumerate(tables.items(), 1):
            S.append(Paragraph(f"Table {ti}. {E(name.replace('_', ' '))}", st_h2))
            cols = list(df.columns)
            max_cols = 9
            if len(cols) <= max_cols:
                chunks = [cols]
            else:
                chunks = [[cols[0]] + cols[i:i + max_cols - 1] for i in range(1, len(cols), max_cols - 1)]
            for ch in chunks:
                data = [[Paragraph(E(c), st_head) for c in ch]]
                for _, row in df[ch].iterrows():
                    data.append([Paragraph(E(_fmt(row[c])[:140]), st_cell) for c in ch])
                weights = []
                for c in ch:
                    longest = max([len(str(c)) * 0.8] + [min(len(_fmt(v)), 40) for v in df[c].head(60)])
                    weights.append(max(6.0, min(longest, 40.0)))
                tot = sum(weights)
                widths = [avail_w * w / tot for w in weights]
                t = Table(data, colWidths=widths, repeatRows=1)
                t.setStyle(TableStyle([
                    ("BACKGROUND", (0, 0), (-1, 0), navy),
                    ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#999999")),
                    ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F3F5F9")]),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 3), ("RIGHTPADDING", (0, 0), (-1, -1), 3),
                ]))
                S += [t, Spacer(1, 0.35 * cm)]

    # ---- limits
    S += [PageBreak(), Paragraph("Limitations and sources", st_h1),
          Paragraph("• Predictions only: no dose, no real cells. Confirm ideas with laboratory work "
                    "(for example caspase or Annexin V assays).", st_body),
          Paragraph("• The antioxidant-to-pro-oxidant switch depends on concentration, metal ions and cell type; "
                    "ProTox-3 gives no dose information.", st_body),
          Paragraph("• Scores in RedoxScope are transparent teaching models, not validated predictors.", st_body),
          Paragraph("• Sources: ProTox-3 (https://tox.charite.de/protox3/) and CB-Dock2 "
                    "(http://cao.labshare.cn:11080/cb-dock2/). Please cite them in any report that uses these results.",
                    st_body)]
    doc.build(S)
    return buf.getvalue()
