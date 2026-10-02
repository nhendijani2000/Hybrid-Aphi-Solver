"""Render a project Markdown document as a styled, self-contained HTML report.

    python md2report.py docs/GAUGE_CHOICE.md docs/GaugeChoice.html

WHY THIS EXISTS. The case reports carry a hand-written page.html beside their
.md, which works because each figure needs its own caption and layout. A prose
document has no such need, and the first edition of docs/GAUGE_CHOICE.md was
maintained as two hand-written files that immediately drifted apart. One source,
one renderer.

Deliberately a small converter, not a Markdown implementation: it handles the
constructs this project's prose documents actually use -- ATX headings, fenced
code, pipe tables, blockquotes, ordered/unordered lists, and inline code, bold,
italic and links. Anything else passes through as a paragraph. If a document
needs more than this, the right fix is usually to simplify the document.

The stylesheet is the same one the case reports use, so a docs/ report and a
regression-test report look like the same publication.
"""
import html
import io
import os
import re
import sys

STYLE = """<style>
  :root{
    --ink:#15181d; --ink-2:#4a5260; --ink-3:#6f7888;
    --bg:#fbfaf8; --card:#ffffff; --rule:#e3e0da;
    --accent:#8a3324; --accent-soft:#f6ece9;
    --ok:#2c5f2d; --warn:#8a5a00;
    color-scheme: light;
  }
  @media (prefers-color-scheme: dark){
    :root:not([data-theme="light"]){
      --ink:#e9eaee; --ink-2:#a8b0bd; --ink-3:#7e8796;
      --bg:#14161a; --card:#1b1e24; --rule:#2c313a;
      --accent:#e08b74; --accent-soft:#2a1f1c;
      --ok:#8fc98f; --warn:#e0b467;
      color-scheme: dark;
    }
  }
  :root[data-theme="dark"]{
    --ink:#e9eaee; --ink-2:#a8b0bd; --ink-3:#7e8796;
    --bg:#14161a; --card:#1b1e24; --rule:#2c313a;
    --accent:#e08b74; --accent-soft:#2a1f1c;
    --ok:#8fc98f; --warn:#e0b467;
    color-scheme: dark;
  }
  body{background:var(--bg); color:var(--ink);
    font-family:"IBM Plex Sans",system-ui,sans-serif; line-height:1.6; margin:0;}
  .wrap{max-width:900px; margin:0 auto; padding-inline:16px; padding-block:40px 72px;}
  h1{font-family:"Source Serif 4",Georgia,serif; font-weight:600;
     font-size:clamp(1.9rem,5vw,2.6rem); line-height:1.15; margin:0 0 .5em;
     text-wrap:balance; letter-spacing:-.01em;}
  h2{font-family:"Source Serif 4",Georgia,serif; font-weight:600; font-size:1.5rem;
     margin:2.6rem 0 .2em; padding-top:1.5rem; border-top:1px solid var(--rule);
     text-wrap:balance;}
  h3{font-size:1.02rem; font-weight:600; margin:1.8rem 0 .3em; color:var(--ink);}
  p{margin:.7em 0; max-width:68ch;}
  code{font-family:"IBM Plex Mono",monospace; font-size:.88em;
       background:var(--accent-soft); color:var(--accent);
       padding:.08em .34em; border-radius:3px;}
  pre{background:var(--card); border:1px solid var(--rule); border-radius:6px;
      padding:.8rem 1rem; overflow-x:auto; font-family:"IBM Plex Mono",monospace;
      font-size:.82rem; line-height:1.5; margin:1rem 0;}
  pre code{background:none; color:var(--ink); padding:0; font-size:1em;}
  .tw{overflow-x:auto; margin:1.1rem 0;}
  table{border-collapse:collapse; width:100%; font-size:.88rem;
        font-variant-numeric:tabular-nums;}
  th,td{text-align:left; padding:.45rem .7rem; border-bottom:1px solid var(--rule);
        vertical-align:top;}
  th{font-weight:600; color:var(--ink-2); font-size:.78rem;
     text-transform:uppercase; letter-spacing:.04em;}
  tbody tr:last-child td{border-bottom:none;}
  blockquote{background:var(--card); border:1px solid var(--rule);
      border-left:3px solid var(--accent); border-radius:6px;
      padding:.9rem 1.1rem; margin:1.3rem 0;}
  blockquote p{margin:.35em 0;}
  blockquote strong{color:var(--accent);}
  ul,ol{max-width:68ch; padding-left:1.15rem;}
  li{margin:.35em 0;}
  hr{border:none; border-top:1px solid var(--rule); margin:2.2rem 0;}
  a{color:var(--accent);}
  @media print{
    :root{ color-scheme: light; }
    body{ background:#fff; -webkit-print-color-adjust:exact;
          print-color-adjust:exact; orphans:3; widows:3; }
    .wrap{ max-width:none; padding-block:0 .5rem; }
    h1{ font-size:2rem; }
    h2{ margin-top:1.2rem; padding-top:.7rem; font-size:1.3rem;
        break-after:avoid; page-break-after:avoid; }
    h3{ margin-top:.9rem; break-after:avoid; page-break-after:avoid; }
    .tw, table{ break-inside:avoid; page-break-inside:avoid; }
    tr{ break-inside:avoid; page-break-inside:avoid; }
    p, li, blockquote, ul, ol, pre{ break-inside:auto; page-break-inside:auto; }
    a{ color:inherit; text-decoration:none; }
  }
</style>"""

HEAD = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>%s</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Source+Serif+4:opsz,wght@8..60,400;8..60,600&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
%s</head><body>
<div class="wrap">
"""


def inline(s):
    """Inline spans. Code first, so its contents are not re-processed."""
    out, parts = [], re.split(r"(`[^`]+`)", s)
    for p in parts:
        if p.startswith("`") and p.endswith("`") and len(p) > 1:
            out.append("<code>%s</code>" % html.escape(p[1:-1]))
            continue
        p = html.escape(p)
        p = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', p)
        p = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", p)
        p = re.sub(r"(?<![*\w])\*([^*]+)\*(?!\w)", r"<em>\1</em>", p)
        out.append(p)
    return "".join(out)


def render_table(rows):
    cells = [[c.strip() for c in r.strip().strip("|").split("|")] for r in rows]
    # Row 1 is the header, row 2 the alignment rule, the rest the body.
    body_from = 2 if len(cells) > 1 and set("".join(cells[1])) <= set("-: ") else 1
    head = cells[0]
    out = ['<div class="tw"><table>', "<thead><tr>"]
    out += ["<th>%s</th>" % inline(c) for c in head]
    out += ["</tr></thead>", "<tbody>"]
    for row in cells[body_from:]:
        out.append("<tr>" + "".join("<td>%s</td>" % inline(c) for c in row) + "</tr>")
    out += ["</tbody></table></div>"]
    return "\n".join(out)


def convert(md):
    lines = md.split("\n")
    out, i, title = [], 0, None
    while i < len(lines):
        ln = lines[i]

        if ln.startswith("```"):                      # fenced code
            i += 1
            buf = []
            while i < len(lines) and not lines[i].startswith("```"):
                buf.append(lines[i])
                i += 1
            i += 1
            out.append("<pre><code>%s</code></pre>" % html.escape("\n".join(buf)))
            continue

        if ln.startswith("|"):                        # pipe table
            buf = []
            while i < len(lines) and lines[i].startswith("|"):
                buf.append(lines[i])
                i += 1
            out.append(render_table(buf))
            continue

        if ln.startswith(">"):                        # blockquote
            buf = []
            while i < len(lines) and (lines[i].startswith(">") or
                                      (buf and lines[i].strip() == "")):
                if lines[i].strip() == "":
                    if i + 1 < len(lines) and lines[i + 1].startswith(">"):
                        buf.append("")
                        i += 1
                        continue
                    break
                buf.append(re.sub(r"^>\s?", "", lines[i]))
                i += 1
            inner, para = [], []
            for b in buf + [""]:
                if b.strip() == "":
                    if para:
                        inner.append("<p>%s</p>" % inline(" ".join(para)))
                        para = []
                else:
                    para.append(b.strip())
            out.append("<blockquote>%s</blockquote>" % "".join(inner))
            continue

        m = re.match(r"^(#{1,4})\s+(.*)$", ln)        # heading
        if m:
            lvl, txt = len(m.group(1)), m.group(2).strip()
            if lvl == 1 and title is None:
                title = re.sub(r"<[^>]+>", "", txt)
            out.append("<h%d>%s</h%d>" % (lvl, inline(txt), lvl))
            i += 1
            continue

        if re.match(r"^(---+|\*\*\*+)\s*$", ln):      # rule
            out.append("<hr>")
            i += 1
            continue

        m = re.match(r"^(\s*)([-*]|\d+\.)\s+(.*)$", ln)   # list
        if m:
            ordered = bool(re.match(r"\d+\.", m.group(2)))
            tag = "ol" if ordered else "ul"
            items = []
            while i < len(lines):
                m2 = re.match(r"^(\s*)([-*]|\d+\.)\s+(.*)$", lines[i])
                if m2:
                    items.append(m2.group(3).strip())
                    i += 1
                elif lines[i].strip() and lines[i].startswith(("   ", "\t")) and items:
                    items[-1] += " " + lines[i].strip()     # continuation
                    i += 1
                else:
                    break
            out.append("<%s>%s</%s>"
                       % (tag, "".join("<li>%s</li>" % inline(x) for x in items), tag))
            continue

        if ln.strip() == "":
            i += 1
            continue

        buf = []                                      # paragraph
        while i < len(lines) and lines[i].strip() and not re.match(
                r"^(#{1,4}\s|\||>|```|---+\s*$|\s*([-*]|\d+\.)\s)", lines[i]):
            buf.append(lines[i].strip())
            i += 1
        out.append("<p>%s</p>" % inline(" ".join(buf)))

    return title or "Report", "\n".join(out)


if __name__ == "__main__":
    src, dst = sys.argv[1], sys.argv[2]
    title, body = convert(io.open(src, encoding="utf-8").read())
    io.open(dst, "w", encoding="utf-8", newline="\n").write(
        (HEAD % (html.escape(title), STYLE)) + body + "\n</div>\n</body></html>\n")
    print("  %s -> %s  (%d bytes, title %r)"
          % (os.path.basename(src), os.path.basename(dst),
             os.path.getsize(dst), title))
