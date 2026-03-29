# TailfindR compatibility

> **Note:** TailfindR compatibility applies only to the **Guppy legacy
> pipeline**
> ([`check_tails_guppy()`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)).

Since version 1.0.2, **ninetails** is compatible with
[tailfindR](https://github.com/adnaniazi/tailfindr).

## Recommended approach

We still recommend `nanopolish polya` for determining tail coordinates
because:

- The ninetails model was trained on Nanopolish output
- Nanopolish provides quality metrics that ninetails uses for filtering
- Nanopolish manages resources better and works faster

However, for users who prefer tailfindR, we provide compatibility
functions.

## Using ninetails with tailfindR output

### Step 1: Convert tailfindR output

Before running the pipeline, convert tailfindR output to
ninetails-compatible format:

``` r
converted_tailfindr <- ninetails::convert_tailfindr_output(
  tailfindr_output = '/path/to/tailfindr/output.csv'
)
```

### Step 2: Run the pipeline

Pass the converted output to
[`check_tails_guppy()`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
as the `polya_data` argument:

``` r
results <- ninetails::check_tails_guppy(
  polya_data = converted_tailfindr,
  sequencing_summary = '/path/to/sequencing_summary.txt',
  workspace = '/path/to/workspace',
  num_cores = 2,
  basecall_group = 'Basecall_1D_000',
  pass_only = TRUE,
  save_dir = '~/output/'
)
```

## Important caveats

### No quality metrics

TailfindR does not provide signal quality metrics like Nanopolish does.
This means:

- Signals with poor quality may be included in the analysis
- With Nanopolish, such signals would be filtered based on the quality
  tag
- Results from tailfindR should be treated with additional caution

### Not seamlessly integrated

Due to the quality metric limitation, ninetails has not been seamlessly
integrated with tailfindR like it has with Nanopolish. The explicit
conversion step is required.

## Comparison

| Feature               | Nanopolish    | TailfindR           |
|-----------------------|---------------|---------------------|
| Quality metrics       | ✅ Yes        | ❌ No               |
| Resource management   | Better        | Variable            |
| Speed                 | Faster        | Variable            |
| Ninetails integration | Seamless      | Requires conversion |
| Recommendation        | **Preferred** | Use with caution    |

## Summary

For best results with the Guppy legacy pipeline:

1.  **Preferred**: Use `nanopolish polya` output directly
2.  **Alternative**: Convert tailfindR output using
    [`convert_tailfindr_output()`](https://LRB-IIMCB.github.io/ninetails/reference/convert_tailfindr_output.md),
    but interpret results cautiously

For new analyses, we recommend using the **Dorado DRS pipeline**
([`check_tails_dorado_DRS()`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_dorado_DRS.md))
instead of the Guppy legacy pipeline.
