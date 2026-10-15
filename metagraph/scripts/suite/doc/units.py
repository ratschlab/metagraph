"""Units: the reading model of the results document (DESIGN-DOC §2, §10). A scenario belongs to exactly one unit,
a unit to one family, a family to one chapter. The metadata comes from the manifest first (a `doc` block per
scenario: v5) and else from a mapping file (`doc/units/suite_v<version>.json`: v4). `resolve` builds the model,
`validate` lists what is wrong with it (the maintenance test fails on any problem), `prose_scaffold` writes the
placeholder entries a prose file lacks.
"""

import json
import os

UNIT_TYPES = ("series", "replicates", "refusals", "list", "pairs")
CHAPTER_ORDER = ["R", "T", "P", "L", "X", "M", "Q", "C", "N", "D"]
CHAPTER_TITLES = {"R": "Resolve", "T": "Traverse", "P": "Patterns up to k", "L": "Patterns longer than k",
                  "X": "Peptides", "M": "Mismatches", "Q": "Predicates and motif",
                  "C": "Capabilities and protocol", "N": "Negatives and near-misses", "D": "Determinism"}
SERIES_MAX = 8
HERE = os.path.dirname(os.path.abspath(__file__))


class UnitsError(Exception):
    """Unit metadata the generator cannot use (exit 2)."""


def default_mapping_path(manifest):
    version = (manifest or {}).get("manifest_version")
    return os.path.join(HERE, "units", "suite_v%s.json" % version)


def load_mapping(path):
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict) or "_units" not in data:
        raise UnitsError("%s: not a unit mapping (no _units)" % path)
    return data


class Units(object):
    """chapters: [letter] in document order; titles: {letter: title}; families: {fid: {chapter, title, units}};
    units: {uid: {family, unit_type, unit_title, property, members, groups, relations}}; by_scenario: {sid: {chapter,
    family, unit_id, unit_type, variant, group}}; relation_unit: {relation id: unit id}."""

    def __init__(self):
        self.chapters = []
        self.titles = {}
        self.families = {}
        self.units = {}
        self.by_scenario = {}
        self.relation_unit = {}
        self.source = None

    def chapter_families(self, chapter):
        return [fid for fid, f in self.families.items() if f.get("chapter") == chapter]

    def chapter_scenarios(self, chapter):
        return [sid for sid, m in self.by_scenario.items() if m.get("chapter") == chapter]

    def family_scenarios(self, fid):
        out = []
        for uid in (self.families.get(fid) or {}).get("units") or []:
            out.extend((self.units.get(uid) or {}).get("members") or [])
        return out

    def unit_of(self, sid):
        return (self.by_scenario.get(sid) or {}).get("unit_id")


def _chapter_of(sid, meta):
    return meta.get("chapter") or (sid.split("-")[0] if "-" in sid else sid[:1])


def resolve(manifest, mapping=None):
    """The model from the manifest's `doc` blocks, else from `mapping`; a scenario in neither raises UnitsError
    naming the ids."""
    units = Units()
    scenarios = [s for s in (manifest or {}).get("scenarios") or [] if isinstance(s, dict) and s.get("id")]
    top_units = (manifest or {}).get("units") if isinstance((manifest or {}).get("units"), dict) else {}
    mapping = mapping or {}
    unmapped = []
    for s in scenarios:
        sid = s["id"]
        doc = s.get("doc") if isinstance(s.get("doc"), dict) else None
        if doc and doc.get("unit"):
            meta = {"chapter": doc.get("chapter") or s.get("area") or _chapter_of(sid, {}),
                    "family": doc.get("family") or sid.rsplit("-", 1)[0], "unit_id": doc["unit"],
                    "unit_type": doc.get("unit_type") or "list", "variant": doc.get("variant") or s.get("title"),
                    "group": doc.get("group"), "unit_title": doc.get("unit_title") or doc["unit"],
                    "relations": doc.get("relations") or [], "family_title": doc.get("family_title")}
            units.source = units.source or "manifest"
        elif isinstance(mapping.get(sid), dict):
            m = mapping[sid]
            meta = {"chapter": m.get("chapter") or s.get("area") or _chapter_of(sid, {}), "family": m.get("family"),
                    "unit_id": m.get("unit_id"), "unit_type": m.get("unit_type") or "list",
                    "variant": m.get("variant") or s.get("title"), "group": m.get("group"),
                    "unit_title": m.get("unit_title"), "relations": [], "family_title": None}
            units.source = units.source or "mapping"
        else:
            unmapped.append(sid)
            continue
        units.by_scenario[sid] = meta
        uid = meta["unit_id"]
        u = units.units.setdefault(uid, {"family": meta["family"], "unit_type": meta["unit_type"],
                                         "unit_title": meta["unit_title"], "property": None, "members": [],
                                         "groups": None, "relations": []})
        u["members"].append(sid)
        for rid in meta.get("relations") or []:
            if rid not in u["relations"]:
                u["relations"].append(rid)
        fid = meta["family"]
        f = units.families.setdefault(fid, {"chapter": meta["chapter"], "title": meta.get("family_title") or fid,
                                            "units": []})
        if uid not in f["units"]:
            f["units"].append(uid)
    if unmapped:
        raise UnitsError("%d scenario(s) without unit metadata (no doc block, not in the mapping): %s" % (
            len(unmapped), ", ".join(unmapped[:20]) + (" ..." if len(unmapped) > 20 else "")))
    # the mapping's own unit and family entries (titles, property lines, groups, relations, unit order)
    for uid, info in (mapping.get("_units") or {}).items():
        if uid in units.units and isinstance(info, dict):
            u = units.units[uid]
            u["unit_title"] = info.get("unit_title") or u["unit_title"]
            u["property"] = info.get("property")
            u["groups"] = info.get("groups")
            for rid in info.get("relations") or []:
                if rid not in u["relations"]:
                    u["relations"].append(rid)
    for fid, info in (mapping.get("_families") or {}).items():
        if fid in units.families and isinstance(info, dict):
            f = units.families[fid]
            f["title"] = info.get("title") or f["title"]
            ordered = [x for x in (info.get("units") or []) if x in f["units"]]
            f["units"] = ordered + [x for x in f["units"] if x not in ordered]
    for uid, info in top_units.items():
        if uid in units.units and isinstance(info, dict):
            u = units.units[uid]
            for rid in info.get("relations") or []:
                if rid not in u["relations"]:
                    u["relations"].append(rid)
            if info.get("groups_order"):
                u["groups"] = info["groups_order"]
    for u in units.units.values():
        for rid in u["relations"]:
            units.relation_unit.setdefault(rid, next(k for k, v in units.units.items() if v is u))
    # chapters in the fixed order, then any other letter in order of appearance
    seen = []
    for sid in units.by_scenario:
        ch = units.by_scenario[sid]["chapter"]
        if ch not in seen:
            seen.append(ch)
    units.chapters = [c for c in CHAPTER_ORDER if c in seen] + [c for c in seen if c not in CHAPTER_ORDER]
    titles = dict(CHAPTER_TITLES)
    titles.update((mapping.get("_chapters") or {}) if isinstance(mapping.get("_chapters"), dict) else {})
    units.titles = dict((c, titles.get(c, c)) for c in units.chapters)
    return units


def validate(units, manifest):
    """The problems of a model: every scenario once, units non-empty, types in the vocabulary, series of at most
    SERIES_MAX distinct variants, groups consistent, every relation under exactly one unit whose members include a
    or b, every family's units and chapter present."""
    problems = []
    scenarios = [s["id"] for s in (manifest or {}).get("scenarios") or [] if isinstance(s, dict) and s.get("id")]
    counted = {}
    for uid, u in units.units.items():
        if not u["members"]:
            problems.append("unit %s has no member" % uid)
        if u["unit_type"] not in UNIT_TYPES:
            problems.append("unit %s: type %r not in %s" % (uid, u["unit_type"], ", ".join(UNIT_TYPES)))
        for sid in u["members"]:
            counted[sid] = counted.get(sid, 0) + 1
        if u["unit_type"] == "series":
            variants = [units.by_scenario[s]["variant"] for s in u["members"] if s in units.by_scenario]
            if len(variants) > SERIES_MAX:
                problems.append("series %s has %d variants (at most %d)" % (uid, len(variants), SERIES_MAX))
            if len(set(variants)) != len(variants):
                problems.append("series %s repeats a variant label" % uid)
        if u.get("groups"):
            for sid in u["members"]:
                g = units.by_scenario.get(sid, {}).get("group")
                if g not in u["groups"]:
                    problems.append("unit %s: member %s in group %r not in the unit's groups" % (uid, sid, g))
        if u["family"] not in units.families:
            problems.append("unit %s: family %s unknown" % (uid, u["family"]))
        elif uid not in units.families[u["family"]]["units"]:
            problems.append("family %s does not list its unit %s" % (u["family"], uid))
    for sid in scenarios:
        n = counted.get(sid, 0)
        if n != 1:
            problems.append("scenario %s in %d unit(s)" % (sid, n))
    for fid, f in units.families.items():
        if f.get("chapter") not in units.chapters:
            problems.append("family %s: chapter %r unknown" % (fid, f.get("chapter")))
    members_of = dict((uid, set(u["members"])) for uid, u in units.units.items())
    for rel in (manifest or {}).get("relations") or []:
        if not isinstance(rel, dict) or not rel.get("id"):
            continue
        holders = [uid for uid, u in units.units.items() if rel["id"] in u["relations"]]
        if len(holders) != 1:
            problems.append("relation %s under %d unit(s)" % (rel["id"], len(holders)))
        for uid in holders:
            if rel.get("a") not in members_of[uid] and rel.get("b") not in members_of[uid]:
                problems.append("relation %s under %s, whose members hold neither %s nor %s" % (
                    rel["id"], uid, rel.get("a"), rel.get("b")))
    return problems


def ids_text(members, groups=None, by_group=None):
    """The id line under a table: the common prefix once, then the numbers, ranges where consecutive; by group
    when `by_group` gives {group: [ids]} in group order."""
    import re
    from latex import esc

    def split(sid):
        m = re.match(r"^(.*-)(\d+)$", sid)
        return (m.group(1), m.group(2)) if m else (None, sid)

    def ranges(ids):
        out, run = [], []
        for sid in ids:
            prefix, num = split(sid)
            if run and run[-1][0] == prefix and prefix is not None and int(num) == int(run[-1][1]) + 1 \
                    and len(num) == len(run[-1][1]):
                run.append((prefix, num))
            else:
                if run:
                    out.append(run)
                run = [(prefix, num)]
        if run:
            out.append(run)
        parts = []
        for r in out:
            if len(r) >= 3:
                parts.append("%s–%s" % (r[0][1], r[-1][1]))
            else:
                parts.extend(x[1] for x in r)
        return parts

    prefixes = set(split(s)[0] for s in members)
    if len(prefixes) == 1 and None not in prefixes:
        prefix = prefixes.pop()
        head = "ids: \\texttt{%s}" % esc(prefix)
        if by_group:
            body = " · ".join("%s (%s)" % (" ".join(ranges(ids)), esc(g)) for g, ids in by_group.items() if ids)
        else:
            body = " · ".join(ranges(members))
        return head + body
    return "ids: " + " · ".join("\\texttt{%s}" % esc(s) for s in members)


def prose_scaffold(units, prose, chapter, check_kinds=None):
    """The prose dict of `chapter` with every missing key filled by a placeholder (existing text untouched):
    {chapter, title, intro, families{id: {title, prose}}, units{id: text}, checks{name: text}}."""
    out = dict(prose or {})
    out.setdefault("chapter", chapter)
    out.setdefault("title", units.titles.get(chapter, chapter))
    if not out.get("intro"):
        out["intro"] = ["(chapter prose: to be written)"]
    fams = dict(out.get("families") or {})
    for fid in units.chapter_families(chapter):
        entry = dict(fams.get(fid) or {})
        entry.setdefault("title", units.families[fid].get("title") or fid)
        if not entry.get("prose"):
            entry["prose"] = "(family prose: to be written)"
        fams[fid] = entry
    out["families"] = fams
    us = dict(out.get("units") or {})
    for fid in units.chapter_families(chapter):
        for uid in units.families[fid]["units"]:
            if not us.get(uid):
                us[uid] = units.units[uid].get("property") or "(unit property: to be written)"
    out["units"] = us
    checks = dict(out.get("checks") or {})
    for name in check_kinds or []:
        if not checks.get(name):
            checks[name] = "(check prose: to be written)"
    if checks or "checks" in out:
        out["checks"] = checks
    return out
