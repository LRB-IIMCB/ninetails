# Build ggplot2 annotation layers for non-A residue highlights

Creates a list of
[`ggplot2::annotate()`](https://ggplot2.tidyverse.org/reference/annotate.html)
layers (semi-transparent rectangles + residue labels) for each non-A
modification in a single read. The rectangles are drawn behind the
signal line to highlight the estimated modification position with a
flank.

## Usage

``` r
.build_nonA_overlay(
  read_nonA_data,
  poly_tail_start,
  poly_tail_end,
  nonA_flank = 250,
  x_scale_factor = 1,
  alpha = 0.15
)
```

## Arguments

- read_nonA_data:

  Data frame. Rows from the `nonadenosine_residues` table for a single
  read. Required columns: `prediction` (C/G/U), `est_nonA_pos`,
  `polya_length`.

- poly_tail_start:

  Numeric. Poly(A) start coordinate in raw signal.

- poly_tail_end:

  Numeric. Poly(A) end coordinate in raw signal.

- nonA_flank:

  Numeric. Number of raw signal positions to highlight on each side of
  the estimated modification center. Default 250 (matching the
  ~500-point raw signal extent of a 100-point interpolated chunk).

- x_scale_factor:

  Numeric. Scaling factor for x-axis coordinates. Use 1 for raw position
  plots, `1 / sampling_rate` for time-rescaled plots. Default 1.

- alpha:

  Numeric. Transparency of the highlight rectangles \[0, 1\]. Default
  0.15.

## Value

A list of ggplot2 annotation layers (can be added to a ggplot with `+`).
Returns an empty list if `read_nonA_data` has zero rows or is NULL.

## Details

Residue colors match the package defaults used in
[`plot_residue_counts`](https://LRB-IIMCB.github.io/ninetails/reference/plot_residue_counts.md)
and
[`plot_panel_characteristics`](https://LRB-IIMCB.github.io/ninetails/reference/plot_panel_characteristics.md):
C = `"#3a424f"`, G = `"#50a675"`, U = `"#b0bdd4"`.
