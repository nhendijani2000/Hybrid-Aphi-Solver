# -*- coding: utf-8 -*-
"""Page-accurate PDF text extraction for LaTeX documents with CID fonts.

    pvpython pdftext2.py in.pdf out.txt [first_page last_page]

WHY ex2.py IS NOT ENOUGH. ex2.py only understands literal strings -- `(text) Tj`
and `[(a) -20 (b)] TJ`. A pdftex document with subsetted Type1/TrueType fonts
writes HEX strings of two-byte glyph ids instead:

    [<00340055>-20<0042>18<0043>] TJ        "Stab..."

ex2.py's string regex matches none of that, so Stysch's 120-page thesis yielded
14 k characters from 148 content streams -- roughly the handful of pages that
happened to use literal strings.

Glyph ids are font-internal, NOT character codes. For this thesis the main text
font happens to satisfy `char = glyph + 31`, but that offset is an accident of
one subset's glyph order and is wrong for the math and italic fonts. So instead
of guessing, this reads the /ToUnicode CMap each font carries (195 of them here)
and decodes through it.

It also walks the page tree, so output is split by PAGE with the page's own
index -- necessary to answer "read chapter 6, pages 63-86" rather than hunting
through one undifferentiated blob.

Deliberately still a reader, not a PDF library: no encryption, no shading, no
reconstruction of reading order for multi-column layouts. Text comes out in
content-stream order, which for LaTeX output is reading order.
"""
import io
import re
import sys
import zlib

WS = rb"[\x00\t\n\f\r ]"


# ------------------------------------------------------------- objects ----
def _inflate(raw, dct):
    """Best-effort decode of a stream body."""
    if b"FlateDecode" in dct:
        try:
            return zlib.decompress(raw)
        except zlib.error:
            # Trailing junk or a truncated final block: salvage what we can.
            try:
                return zlib.decompressobj().decompress(raw)
            except zlib.error:
                return None
    return raw


def parse_objects(data):
    """num -> (dict_bytes, stream_bytes|None), including objects in ObjStms."""
    objs = {}
    for m in re.finditer(rb"(\d+)%s+0%s+obj\b" % (WS, WS), data):
        num, start = int(m.group(1)), m.end()
        # The dictionary is small; a stream keyword inside the next 4 kB means
        # this is a stream object. Searching for `endobj` first is unsafe --
        # binary stream data contains that byte sequence often enough.
        window = data[start:start + 4096]
        sm = re.search(rb"stream\r?\n", window)
        if sm:
            dct = window[:sm.start()]
            s0 = start + sm.end()
            s1 = data.find(b"endstream", s0)
            objs[num] = (dct, data[s0:s1] if s1 > 0 else None)
        else:
            e = data.find(b"endobj", start)
            objs[num] = (data[start:e if e > 0 else start + 4096], None)

    # PDF >= 1.5 hides most objects inside compressed object streams.
    for num, (dct, raw) in list(objs.items()):
        if b"/ObjStm" not in dct or raw is None:
            continue
        body = _inflate(raw, dct)
        if body is None:
            continue
        n = re.search(rb"/N%s+(\d+)" % WS, dct)
        first = re.search(rb"/First%s+(\d+)" % WS, dct)
        if not (n and first):
            continue
        n, first = int(n.group(1)), int(first.group(1))
        nums = [int(x) for x in body[:first].split()]
        for i in range(n):
            onum, off = nums[2 * i], nums[2 * i + 1]
            end = nums[2 * i + 3] + first if i + 1 < n else len(body)
            if onum not in objs:                  # top-level wins
                objs[onum] = (body[first + off:end], None)
    return objs


def deref(objs, blob, key):
    """Value of `key` in `blob`, following one level of indirection."""
    m = re.search(key + rb"%s*(\d+)%s+0%s+R" % (WS, WS, WS), blob)
    if m:
        num = int(m.group(1))
        return objs.get(num, (b"", None))[0]
    m = re.search(key + rb"%s*(<<.*?>>|\[.*?\])" % WS, blob, re.S)
    return m.group(1) if m else None


# ------------------------------------------------------------- ToUnicode ----
def parse_cmap(text):
    """CMap stream -> {glyph_code: unicode_string}."""
    out = {}

    def utf16(h):
        try:
            b = bytes.fromhex(h.decode("latin-1").strip())
        except (ValueError, AttributeError):
            return ""
        try:
            return b.decode("utf-16-be", "ignore")
        except Exception:
            return ""

    for blk in re.findall(rb"beginbfchar(.*?)endbfchar", text, re.S):
        for src, dst in re.findall(rb"<([0-9A-Fa-f]+)>%s*<([0-9A-Fa-f]*)>" % WS,
                                   blk):
            out[int(src, 16)] = utf16(dst)

    for blk in re.findall(rb"beginbfrange(.*?)endbfrange", text, re.S):
        # <lo> <hi> <dst>
        for lo, hi, dst in re.findall(
                rb"<([0-9A-Fa-f]+)>%s*<([0-9A-Fa-f]+)>%s*<([0-9A-Fa-f]*)>"
                % (WS, WS), blk):
            lo, hi = int(lo, 16), int(hi, 16)
            base = utf16(dst)
            for k in range(lo, min(hi, lo + 65535) + 1):
                if base:
                    out[k] = base[:-1] + chr(ord(base[-1]) + k - lo)
        # <lo> <hi> [<d1> <d2> ...]
        for lo, _hi, arr in re.findall(
                rb"<([0-9A-Fa-f]+)>%s*<([0-9A-Fa-f]+)>%s*\[(.*?)\]"
                % (WS, WS), blk, re.S):
            lo = int(lo, 16)
            for i, dst in enumerate(re.findall(rb"<([0-9A-Fa-f]*)>", arr)):
                out[lo + i] = utf16(dst)
    return out


def font_maps(objs):
    """font object number -> {glyph: unicode}, for every font with a ToUnicode."""
    maps = {}
    for num, (dct, _) in objs.items():
        if b"/Font" not in dct and b"/BaseFont" not in dct:
            continue
        m = re.search(rb"/ToUnicode%s*(\d+)%s+0%s+R" % (WS, WS, WS), dct)
        if not m:
            continue
        tnum = int(m.group(1))
        tdct, traw = objs.get(tnum, (b"", None))
        if traw is None:
            continue
        body = _inflate(traw, tdct)
        if body:
            maps[num] = parse_cmap(body)
    return maps


# ------------------------------------------------------------- pages ----
def page_objects(data, objs):
    """/Type /Page objects in file order -- document order for pdftex output."""
    pages = []
    for m in re.finditer(rb"(\d+)%s+0%s+obj\b" % (WS, WS), data):
        num = int(m.group(1))
        dct = objs.get(num, (b"", None))[0]
        if re.search(rb"/Type%s*/Page[^s]" % WS, dct):
            pages.append((m.start(), num))
    pages.sort()
    seen, out = set(), []
    for _, num in pages:
        if num not in seen:
            seen.add(num)
            out.append(num)
    if out:
        return out
    # Fallback: objects inside ObjStms have no file offset of their own.
    return [n for n, (d, _) in sorted(objs.items())
            if re.search(rb"/Type%s*/Page[^s]" % WS, d)]


def page_fonts(objs, pdct):
    """resource name (b'F25') -> glyph map, for one page."""
    res = deref(objs, pdct, rb"/Resources")
    if not res:
        return {}
    fdict = deref(objs, res, rb"/Font")
    if not fdict:
        return {}
    out = {}
    for name, onum in re.findall(rb"/([A-Za-z0-9+_.-]+)%s*(\d+)%s+0%s+R"
                                 % (WS, WS, WS), fdict):
        out[name] = int(onum)
    return out


def page_content(objs, pdct):
    """Concatenated, decompressed content streams of one page."""
    m = re.search(rb"/Contents%s*(\d+)%s+0%s+R" % (WS, WS, WS), pdct)
    nums = [int(m.group(1))] if m else []
    if not nums:
        m = re.search(rb"/Contents%s*\[(.*?)\]" % WS, pdct, re.S)
        if m:
            nums = [int(x) for x in
                    re.findall(rb"(\d+)%s+0%s+R" % (WS, WS), m.group(1))]
    parts = []
    for n in nums:
        dct, raw = objs.get(n, (b"", None))
        if raw is not None:
            body = _inflate(raw, dct)
            if body:
                parts.append(body)
    return b"\n".join(parts)


# ------------------------------------------------------------- decode ----
TOKEN = re.compile(
    rb"/([A-Za-z0-9+_.-]+)%s+[\d.]+%s+Tf"      # 1: font switch
    rb"|\[(.*?)\]%s*TJ"                        # 2: array show
    rb"|(<[0-9A-Fa-f\s]*>|\((?:\\.|[^\\()])*\))%s*Tj"   # 3: single show
    rb"|(T\*|Td|TD|TL)"                        # 4: line movers
    rb"|(ET)" % (WS, WS, WS, WS), re.S)

PIECE = re.compile(rb"<[0-9A-Fa-f\s]*>|\((?:\\.|[^\\()])*\)|(-?[\d.]+)", re.S)

# A TJ array interleaves strings with kern adjustments in 1/1000 em. LaTeX emits
# no space glyph: an inter-word gap IS a large negative adjustment. Below this
# threshold the gap is letter kerning, above it a word break. -150 sits well
# clear of both (kerns here run to about -45, word gaps from about -250).
SPACE_KERN = -150.0
ESC = {b"n": b"\n", b"r": b"\r", b"t": b"\t", b"b": b"\b", b"f": b"\f",
       b"(": b"(", b")": b")", b"\\": b"\\"}


def unescape(lit):
    out, i = bytearray(), 1
    while i < len(lit) - 1:
        c = lit[i:i + 1]
        if c == b"\\":
            nxt = lit[i + 1:i + 2]
            if nxt.isdigit():
                j = i + 1
                while j < len(lit) - 1 and lit[j:j + 1].isdigit() and j < i + 4:
                    j += 1
                out.append(int(lit[i + 1:j], 8) & 0xFF)
                i = j
                continue
            out += ESC.get(nxt, nxt)
            i += 2
            continue
        out += c
        i += 1
    return bytes(out)


def show(piece, cmap, two_byte):
    """One string operand -> text."""
    if piece.startswith(b"<"):
        h = re.sub(rb"\s", b"", piece[1:-1])
        if len(h) % 2:
            h += b"0"
        try:
            raw = bytes.fromhex(h.decode("latin-1"))
        except ValueError:
            return ""
        step = 2 if two_byte else 1
        out = []
        for i in range(0, len(raw) - step + 1, step):
            code = int.from_bytes(raw[i:i + step], "big")
            if cmap and code in cmap:
                out.append(cmap[code])
            elif not cmap:
                out.append(chr(code) if 32 <= code < 127 else "")
            else:
                out.append("")
        return "".join(out)
    raw = unescape(piece)
    if cmap:
        return "".join(cmap.get(b, "") for b in raw)
    return raw.decode("latin-1")


def decode_page(content, name2obj, maps):
    cmap, two_byte, out = None, True, []
    for m in TOKEN.finditer(content):
        if m.group(1) is not None:                       # Tf
            onum = name2obj.get(m.group(1))
            cmap = maps.get(onum)
            # A CID font's codes exceed one byte; a simple font's do not.
            two_byte = bool(cmap) and max(cmap) > 0xFF
        elif m.group(2) is not None:                     # TJ
            for p in PIECE.finditer(m.group(2)):
                if p.group(1) is not None:               # a kern adjustment
                    try:
                        if float(p.group(1)) <= SPACE_KERN:
                            out.append(" ")
                    except ValueError:
                        pass
                    continue
                out.append(show(p.group(), cmap, two_byte))
        elif m.group(3) is not None:                     # Tj
            out.append(show(m.group(3), cmap, two_byte))
        elif m.group(4) is not None:                     # new line
            out.append("\n")
        else:                                            # ET
            out.append("\n")
    return re.sub(r"\n{3,}", "\n\n", "".join(out))


def main():
    data = open(sys.argv[1], "rb").read()
    objs = parse_objects(data)
    maps = font_maps(objs)
    pages = page_objects(data, objs)
    lo = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    hi = int(sys.argv[4]) if len(sys.argv) > 4 else len(pages)

    chunks, nchars = [], 0
    for idx, pnum in enumerate(pages, 1):
        if not (lo <= idx <= hi):
            continue
        pdct = objs[pnum][0]
        txt = decode_page(page_content(objs, pdct), page_fonts(objs, pdct), maps)
        nchars += len(txt)
        chunks.append(u"\n\n========== PDF PAGE %d ==========\n%s" % (idx, txt))

    io.open(sys.argv[2], "w", encoding="utf-8").write(u"".join(chunks))
    sys.stdout.write(
        "  %d objects, %d fonts with ToUnicode, %d pages; "
        "wrote pages %d-%d, %d chars\n"
        % (len(objs), len(maps), len(pages), lo, min(hi, len(pages)), nchars))

    # Honesty check: the page tree declares its own count. If we found fewer
    # pages than that, some /Page objects are hiding somewhere this reader does
    # not look (a deeply nested page tree, or an object stream whose container
    # we failed to inflate). Silently returning two thirds of a paper is the
    # one failure mode that would not look like a failure.
    declared = max([int(m.group(1)) for m in
                    re.finditer(rb"/Count%s+(\d+)" % WS, data)] or [0])
    if declared > len(pages):
        sys.stdout.write(
            "  WARNING: page tree declares %d pages, only %d found -- "
            "output is incomplete\n" % (declared, len(pages)))


if __name__ == "__main__":
    main()
