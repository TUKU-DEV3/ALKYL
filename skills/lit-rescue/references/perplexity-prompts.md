# Lit-Rescue — Perplexity Prompt Templates

Templates for querying Perplexity MCP when standard knowledge fails.
Each template is designed to elicit **sourced, verifiable answers** and suppress hallucinated plausibility.

---

## Core Anti-Hallucination Frame

Always wrap any Perplexity query with this frame when the question is technical:

```
[QUERY]

Requirements for your answer:
- Cite only peer-reviewed papers, official documentation, or GitHub issues/changelogs.
- Include DOIs or direct URLs to documentation pages.
- If you are uncertain, say so explicitly — do not guess.
- Do not extrapolate beyond what your sources directly state.
```

---

## Template by Problem Type

### TYPE: METHOD — Which algorithm/approach to use

**Use when:** "How do I compute X?", "Which method for Y?", "Best approach to Z?"

```
I need to [TASK DESCRIPTION] using [LIBRARY/FRAMEWORK if known].

What peer-reviewed methods or established algorithms exist for this?
List approaches in order of typical adoption, with:
- Brief description of each approach
- Known limitations
- DOI or documentation reference for each

Context: [domain, scale of problem, known constraints]

Do not suggest approaches that are not documented in the literature or official docs.
```

**Example expansion:**
```
I need to compute electron density topology (critical points) from a DFT wavefunction using Python.

What peer-reviewed methods or established algorithms exist for this?
List approaches in order of typical adoption, with:
- Brief description of each approach
- Known limitations
- DOI or documentation reference for each

Context: molecular system, single-point DFT output (ORCA), looking for QTAIM analysis.

Do not suggest approaches that are not documented in the literature or official docs.
```

---

### TYPE: PARAM — Parameter values / hyperparameter settings

**Use when:** "What threshold?", "What learning rate?", "Typical cutoff for?"

```
What are the standard / recommended parameter values for [PARAMETER] when doing [TASK]?

I need:
- Typical literature values (with citation)
- How these values vary with [scale / domain / use case]
- Primary reference that established these defaults (DOI preferred)

Domain: [specific field]
Software: [tool + version if relevant]

Only report values that appear in published papers or official documentation.
```

---

### TYPE: BUG — Unexpected library behavior or error

**Use when:** An error occurs that shouldn't based on docs, or behavior contradicts expectations.

```
I am encountering unexpected behavior in [LIBRARY] version [VERSION]:

Problem: [exact description of what happens vs. what is expected]
Code snippet or error message:
[paste here]

Is this a known bug, a version-specific issue, or expected behavior?

Please check:
- Official GitHub issues for [library] matching this behavior
- Changelog entries for [version range]
- StackOverflow or forum discussions about this specific error

Cite the GitHub issue URL or documentation section if found.
```

---

### TYPE: THEORY — Physical/mathematical basis

**Use when:** "Why does X work?", "What is the theoretical justification for Y?"

```
What is the theoretical basis for [CONCEPT/METHOD]?

I need a rigorous explanation, not a simplified analogy.
Please cite the original paper(s) that introduced or formalized this concept (DOI preferred).

Specific question: [formulate precisely]
Domain: [physics / chemistry / ML / biology / etc.]

If there are multiple competing theoretical frameworks, list them and cite each.
```

---

### TYPE: PROTOCOL — Standard experimental or computational procedure

**Use when:** "What is the standard workflow for?", "Step-by-step procedure for?"

```
What is the standard [computational/experimental] protocol for [TASK]?

I need:
- Step-by-step procedure as described in the primary reference
- Key parameters or settings (with typical values)
- Common pitfalls to avoid
- Primary reference paper or protocol database entry (DOI or URL)

Domain: [field]
Constraints: [hardware, software, available inputs]

Cite the methodology paper or protocol that is most widely used in the community.
```

---

### TYPE: BENCHMARK — Comparison of methods / state of the art

**Use when:** "Which is better for?", "State of the art in?", "Comparison of X vs Y?"

```
What does the literature say about the comparative performance of [METHOD A] vs [METHOD B] for [TASK]?

I need:
- Summary of benchmark studies (with DOI)
- Conditions under which each method performs better
- Most recent systematic comparison (prioritize papers from 2022–2026)

Avoid: opinion-based comparisons, blog posts, or unsourced claims.
Domain: [field]
```

---

### TYPE: DOMAIN — Highly specialized / niche question

**Use when:** The question is very specific to a subfield and training knowledge is likely incomplete.

```
[SPECIFIC QUESTION — formulated as precisely as possible]

This is a specialized question in [SUBFIELD].
Please search for:
- Review articles or textbook chapters covering this topic
- Primary research papers with direct evidence
- Official documentation if this involves a software tool

Provide DOIs for all cited papers. If no peer-reviewed source exists for this specific claim, say so explicitly.

Do not fill gaps in your knowledge with plausible-sounding extrapolation.
```

---

## Iteration Strategy

If the first Perplexity query returns insufficient results:

1. **Broaden terminology** — try synonyms or the concept's historical name
   ```
   "QTAIM" → also try "atoms in molecules" "topological analysis" "Bader analysis"
   ```

2. **Change abstraction level** — too specific → go to the parent concept
   ```
   "RadiusOfGyration for flexible polymers in implicit solvent" → "protein radius of gyration MD simulation"
   ```

3. **Split the question** — decompose into two simpler queries
   ```
   "How to use X for Y in context Z" → (1) "How does X work" + (2) "Application of X to Y"
   ```

4. **Target a specific source type** — add `site:docs.[library].org` or `"Journal of Chemical Theory"`

---

## Perplexity MCP — Available Tools

When Perplexity MCP is configured:

```python
# Check availability
mcp__perplexity__search(query="...", focus="academic")  # academic mode preferred

# Focus modes:
# "academic"  → Semantic Scholar, PubMed, arXiv indexed sources
# "internet"  → full web (for bugs, GitHub, docs)
# "writing"   → not useful for lit-rescue
```

If the tool name differs in your setup, check available MCP tools via `mcp__perplexity__*`.

**Preferred call pattern:**
```python
mcp__perplexity__search(
  query="[YOUR FORMULATED QUERY FROM TEMPLATE ABOVE]",
  focus="academic"   # always use academic for METHOD/PARAM/THEORY/BENCHMARK
)
```
