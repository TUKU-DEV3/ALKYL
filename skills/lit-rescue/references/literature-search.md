# Lit-Rescue — PubMed & bioRxiv Search Strategies

Use these MCPs when Perplexity is unavailable or as complementary sources.
Both are always available in the ALKYL environment.

---

## When to Use Each

| Source | Best For | Coverage |
|--------|----------|----------|
| **bioRxiv/medRxiv** | Cutting-edge methods, ML in biology, comp chem, recent benchmarks | Preprints (not peer-reviewed but fast) |
| **PubMed** | Established methods, clinical protocols, biology, biochemistry | Peer-reviewed (slower but vetted) |

**Rule of thumb:**
- Fast-moving field (comp bio, cheminformatics, deep learning) → bioRxiv first
- Established field (enzymology, crystallography, classical MD) → PubMed first
- For maximum coverage → run both, deduplicate by DOI

---

## bioRxiv MCP — Search Patterns

### Basic search
```python
bioRxiv.search_preprints(
    query="[concept] [method] [domain]",
    limit=10
)
```

### Date-filtered search (recent methods)
```python
bioRxiv.search_preprints(
    query="[query]",
    start_date="2022-01-01",
    end_date="2026-03-04",
    limit=10
)
```

### Category-targeted search
```python
# Subject categories available via: bioRxiv.get_categories()
bioRxiv.search_preprints(
    query="[query]",
    category="bioinformatics",  # or: biochemistry, chemistry, pharmacology-toxicology
    limit=10
)
```

### Fetch full preprint details
```python
bioRxiv.get_preprint(doi="10.1101/XXXX.XX.XX.XXXXXX")
# Returns: title, abstract, authors, date, journal (if published)
```

### Check if preprint was published
```python
bioRxiv.search_published_preprints(query="[title keywords]", limit=5)
# Useful to find the peer-reviewed version of a preprint
```

### Query Construction Tips

| Goal | Query pattern |
|------|--------------|
| Find a specific method | `"[method name]" application [domain]` |
| Find benchmarks | `[task] benchmark comparison evaluation 2023 2024` |
| Find negative results | `[method] limitations failure [domain]` |
| Find protocols | `[procedure] protocol step workflow` |
| Track an author's work | `author:[lastname] [topic]` |

---

## PubMed MCP — Search Patterns

### Basic search
```python
PubMed.search_articles(
    query="[concept] [method] [organism/domain]",
    max_results=10,
    sort="relevance"
)
```

### Recent papers only
```python
PubMed.search_articles(
    query="[query]",
    max_results=10,
    sort="date"
)
```

### MeSH-enhanced query (better precision)
```python
# PubMed MeSH terms greatly improve precision
PubMed.search_articles(
    query="[MeSH term][MeSH] AND [method]",
    max_results=10
)
# Example: "Molecular Docking Simulation[MeSH] AND GNINA"
```

### Get full article metadata
```python
PubMed.get_article_metadata(pmid="XXXXXXXX")
# Returns: title, abstract, authors, journal, DOI, MeSH terms
```

### Get full text (if available open access)
```python
PubMed.get_full_text_article(pmid="XXXXXXXX")
```

### Find related articles
```python
PubMed.find_related_articles(pmid="XXXXXXXX", max_results=5)
# Useful to discover follow-up methods papers
```

### ID conversion (PMID ↔ DOI ↔ PMC)
```python
PubMed.convert_article_ids(ids=["PMID:XXXXXX"], target_format="doi")
```

### Query Construction Tips

| Goal | Query pattern |
|------|--------------|
| Specific method | `"[method name]"[Title/Abstract]` |
| Review articles | `[topic] AND "Review"[Publication Type]` |
| Protocols | `[procedure] AND ("protocol"[Title] OR "methods"[Title])` |
| Benchmarks | `[task] AND (benchmark OR comparison OR evaluation OR "systematic review")` |
| Recent only | `[topic] AND ("2022/01/01"[PDAT] : "2026/12/31"[PDAT])` |
| High-impact only | `[topic] AND "Nature"[Journal] OR "Science"[Journal] OR "Cell"[Journal]` |

---

## Combined Search Strategy

For a thorough lit-rescue on a METHOD or BENCHMARK question:

```
1. bioRxiv.search_preprints(query="[method] [domain] benchmark", date_range="2023-2026", limit=10)
   → Find recent preprints (often ahead of publication by 6-18 months)

2. PubMed.search_articles(query="[method] [domain] evaluation", sort="relevance", max_results=10)
   → Find peer-reviewed validation

3. For each promising hit: get_preprint(doi) or get_article_metadata(pmid)
   → Read abstract, extract key findings

4. Cross-check: bioRxiv.search_published_preprints(query="[method name]")
   → Find if preprints were published (adds peer-review weight)
```

---

## Interpreting Results

### Signals of a strong method source
- Paper is in a methods-focused journal (J. Chem. Theory Comput., Bioinformatics, J. Cheminformatics, JCTC, JCIM)
- Paper has >50 citations within 3 years
- Multiple independent groups validated or used the method
- An open-source implementation is provided

### Signals of caution
- Single-group result with no independent replication
- Preprint >2 years old with no journal publication
- Benchmark is self-reported (authors benchmark their own method)
- Dataset is proprietary or the exact dataset is unavailable for comparison

### Reporting format
```
Source: [Authors, Year, Journal, DOI]
Method described: [1-2 sentence summary]
Confidence: ★★★ / ★★☆ / ★☆☆
Applicable because: [why this source answers the specific question]
Limitation: [any caveats about applicability]
```

---

## Fallback: Lookup by Citation

If the user mentions a specific paper title or author+year:

```python
PubMed.lookup_article_by_citation(
    title="[paper title]",
    first_author="[lastname]",
    year=2023
)
```

If a DOI is known:
```python
PubMed.convert_article_ids(ids=["DOI:10.XXXX/XXXXX"], target_format="pmid")
# Then: get_article_metadata(pmid=...)
```
