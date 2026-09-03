#!/usr/bin/env python3
"""Token accounting for the project's Claude Code transcripts.

Sources (all under ~/.claude/projects/<slug>/):
  * <session>.jsonl                       main interactive sessions
  * <session>/subagents/**.jsonl          sub-agent (Agent tool) transcripts, where persisted
  * <session>/workflows/**.jsonl          multi-agent Workflow transcripts, where persisted
  * plus ~/.claude/tasks/<session>/**     task outputs (scanned for .jsonl as well)
Assistant messages are streamed as several records sharing one message.id; usage is de-duplicated per
message id (last record wins).  A second, independent estimate for sub-agents comes from the
<usage><subagent_tokens>..</subagent_tokens>..</usage> blocks that the harness writes into the MAIN
transcript when a sub-agent finishes (de-duplicated per session on the full block text).
"""
import json, glob, os, re, collections
import sys
HOME = os.path.expanduser("~")
# Claude Code stores transcripts under ~/.claude/projects/<slug>/, where <slug> is the project's absolute
# path with '/' replaced by '-'.  Pass the slug (or the full directory) as the first argument.
if len(sys.argv) < 2:
    sys.exit("usage: analyze_transcripts.py <project-slug or ~/.claude/projects/<slug> directory>")
P = sys.argv[1] if os.path.isdir(sys.argv[1]) else HOME + "/.claude/projects/" + sys.argv[1]
TASKS = HOME + "/.claude/tasks"

def scan(path):
    """return dict: model -> Counter(messages, output, thinking, input, cache_read, cache_create), plus tool count, first/last ts"""
    msgs = {}; tools = 0; first = last = None
    with open(path, errors="replace") as fh:
        for line in fh:
            try: r = json.loads(line)
            except Exception: continue
            ts = r.get("timestamp")
            if ts: first = first or ts; last = ts
            if r.get("type") != "assistant": continue
            m = r.get("message", {})
            if isinstance(m.get("content"), list):
                tools += sum(1 for c in m["content"] if isinstance(c, dict) and c.get("type") == "tool_use")
            u = m.get("usage")
            if u: msgs[m.get("id") or r.get("uuid")] = (m.get("model"), u)
    out = collections.defaultdict(collections.Counter)
    prev = None
    for model, u in msgs.values():
        c = out[str(model)]
        ctx = (u.get("input_tokens", 0) or 0) + (u.get("cache_read_input_tokens", 0) or 0) + (u.get("cache_creation_input_tokens", 0) or 0)
        c["incremental"] += ctx if prev is None else max(0, ctx - prev); prev = ctx
        c["messages"] += 1; c["output"] += u.get("output_tokens", 0) or 0
        c["thinking"] += (u.get("output_tokens_details") or {}).get("thinking_tokens", 0) or 0
        c["input"] += u.get("input_tokens", 0) or 0
        c["cache_read"] += u.get("cache_read_input_tokens", 0) or 0
        c["cache_create"] += u.get("cache_creation_input_tokens", 0) or 0
    return out, tools, first, last

def add(total, part):
    for model, c in part.items(): total[model].update(c)

def qwen(model): return "qwen" in model.lower()

# ---------------------------------------------------------------- main sessions
main_files = sorted(glob.glob(P + "/*.jsonl"))
tot_main = collections.defaultdict(collections.Counter); tools_main = 0
print(f"{'first':20} {'last':20} {'session':8} {'models':34} {'msgs':>5} {'tools':>5} {'output':>9} {'input':>11} {'cache_rd':>11}")
for f in main_files:
    part, tools, first, last = scan(f); add(tot_main, part); tools_main += tools
    if part:
        c = sum(part.values(), collections.Counter())
        print(f"{(first or '')[:19]:20} {(last or '')[:19]:20} {os.path.basename(f)[:8]:8} {','.join(sorted(part))[:34]:34} {c['messages']:5} {tools:5} {c['output']:9,} {c['input']:11,} {c['cache_read']:11,}")

# ---------------------------------------------------------------- sub-agent / workflow transcripts
tot_sub = collections.defaultdict(collections.Counter); n_sub_files = 0; tools_sub = 0
per_session_sub = collections.defaultdict(lambda: collections.Counter())
sub_files = sorted(glob.glob(P + "/*/subagents/**/*.jsonl", recursive=True) + glob.glob(P + "/*/workflows/**/*.jsonl", recursive=True)
                   + glob.glob(TASKS + "/**/*.jsonl", recursive=True))
seen = set()
for f in sub_files:
    rp = os.path.realpath(f)
    if rp in seen: continue
    seen.add(rp)
    part, tools, first, last = scan(f)
    if not part: continue
    n_sub_files += 1; tools_sub += tools; add(tot_sub, part)
    sess = f.split(P + "/")[1].split("/")[0] if f.startswith(P) else f.split(TASKS + "/")[1].split("/")[0]
    c = sum(part.values(), collections.Counter()); c["files"] = 1; c["tools"] = tools
    per_session_sub[sess[:8]].update(c)
print(f"\nSUB-AGENT / WORKFLOW TRANSCRIPTS ON DISK: {n_sub_files} files with usage data (of {len(seen)} scanned)")
print(f"{'session':8} {'files':>5} {'msgs':>6} {'tools':>6} {'output':>10} {'input':>12} {'cache_rd':>12}")
for s, c in sorted(per_session_sub.items()):
    print(f"{s:8} {c['files']:5} {c['messages']:6} {c['tools']:6} {c['output']:10,} {c['input']:12,} {c['cache_read']:12,}")

# ---------------------------------------------------------------- embedded completion blocks in main transcripts
pat = re.compile(r"<usage><subagent_tokens>(\d+)</subagent_tokens><tool_uses>(\d+)</tool_uses><duration_ms>(\d+)</duration_ms></usage>")
print("\nSUB-AGENT COMPLETION BLOCKS EMBEDDED IN MAIN TRANSCRIPTS (independent estimate; each block = one finished sub-agent):")
print(f"{'session':8} {'agents':>6} {'subagent_tokens':>16} {'tool_uses':>9} {'hours':>7}")
emb_total = collections.Counter()
for f in main_files:
    txt = open(f, errors="replace").read()
    blocks = set(pat.findall(txt))
    if not blocks: continue
    tok = sum(int(b[0]) for b in blocks); tu = sum(int(b[1]) for b in blocks); ms = sum(int(b[2]) for b in blocks)
    emb_total.update({"agents": len(blocks), "tokens": tok, "tools": tu, "ms": ms})
    print(f"{os.path.basename(f)[:8]:8} {len(blocks):6} {tok:16,} {tu:9} {ms/3.6e6:7.1f}")
print(f"{'TOTAL':8} {emb_total['agents']:6} {emb_total['tokens']:16,} {emb_total['tools']:9} {emb_total['ms']/3.6e6:7.1f}")
print("  (subagent_tokens as reported by the harness; whether it counts output only or input+output is not documented here)")

# ---------------------------------------------------------------- grand totals
def show(title, tot):
    print(f"\n{title}")
    for model, c in sorted(tot.items(), key=lambda kv: -kv[1]["output"]):
        if c["messages"]:
            print(f"  {model:32} messages={c['messages']:6,}  output={c['output']:11,}  input={c['input']:14,}  cache_read={c['cache_read']:14,}  cache_create={c['cache_create']:12,}")
show("MAIN SESSIONS BY MODEL", tot_main)
show("SUB-AGENT/WORKFLOW TRANSCRIPTS BY MODEL", tot_sub)
q_main = sum((c for m, c in tot_main.items() if qwen(m)), collections.Counter())
q_sub = sum((c for m, c in tot_sub.items() if qwen(m)), collections.Counter())
print("\nQWEN PHASE GRAND TOTAL (main + persisted sub-agent transcripts):")
for k in ("messages", "output", "input", "cache_read", "incremental"):
    print(f"  {k:12} main={q_main[k]:14,}  sub-agents={q_sub[k]:14,}  total={q_main[k]+q_sub[k]:14,}")
print(f"  tool calls   main={tools_main:14,}  sub-agents={tools_sub:14,}  total={tools_main+tools_sub:14,}")
tot_ctx = q_main["input"]+q_main["cache_read"]+q_sub["input"]+q_sub["cache_read"]; inc = q_main["incremental"]+q_sub["incremental"]
print(f"\nPREFILL BOUNDS (Qwen): full context re-read every turn = {tot_ctx:,} tokens; incremental (perfect prefix caching) = {inc:,} tokens ({100*inc/tot_ctx:.1f}%)")
for rate in (1000, 2000, 3000):
    print(f"  at {rate} tok/s prefill: no-cache {tot_ctx/rate/3600:6.1f} h   perfect-cache {inc/rate/3600:6.1f} h")
