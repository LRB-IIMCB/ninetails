# Draws an entire squiggle for given read.

Creates segmented plot of raw/rescaled ONT RNA signal (with or without
moves). A standalone function; does not rely on any other preprocessing,
depends solely on Nanopolish, Guppy and fast5 input.

## Usage

``` r
plot_squiggle_fast5(
  readname,
  nanopolish,
  sequencing_summary,
  workspace,
  basecall_group = "Basecall_1D_000",
  moves = FALSE,
  rescale = FALSE
)
```

## Arguments

- readname:

  character string. Name of the given read within the analyzed dataset.

- nanopolish:

  character string. Full path of the .tsv file produced by nanopolish
  polya function.

- sequencing_summary:

  character string. Full path of the .txt file with sequencing summary.

- workspace:

  character string. Full path of the directory to search the basecalled
  fast5 files in. The Fast5 files have to be multi-fast5 file.

- basecall_group:

  character string ("Basecall_1D_000" is set as a default). Name of the
  level in the Fast5 file hierarchy from which the data should be
  extracted.

- moves:

  logical \[TRUE/FALSE\]. If TRUE, moves would be plotted in the
  background as vertical bars/gaps (for values 0/1, respectively) and
  the signal (squiggle) would be plotted in the foreground. Otherwise,
  only the signal would be plotted. As a default, "FALSE" value is set.

- rescale:

  logical \[TRUE/FALSE\]. If TRUE, the signal will be rescaled for
  picoamps (pA) per second (s). If FALSE, raw signal per position will
  be plotted. As a default, the "FALSE" value is set.

## Value

ggplot2 object with squiggle plot depicting nanopore read signal.

## Details

The output plot includes an entire squiggle corresponding to the given
ONT read. Vertical lines mark the 5' (navy blue) and 3' (red) termini of
polyA tail according to the Nanopolish polyA function. In order to
maintain readability of the graph (and to avoid plotting high cliffs -
e.g. jets of the signal caused by a sudden surge of current in the
sensor) the signal is winsorised.

Moves may be plotted only for reads basecalled by Guppy basecaller.
Otherwise the function will throw an error.

## Examples

``` r
if (FALSE) { # \dontrun{

plot <- ninetails::plot_squiggle_fast5(
 readname = "0226b5df-f9e5-4774-bbee-7719676f2ceb",
 nanopolish = system.file('extdata',
                          'test_data',
                          'nanopolish_output.tsv',
                          package = 'ninetails'),
 sequencing_summary = system.file('extdata',
                                  'test_data',
                                  'sequencing_summary.txt',
                                  package = 'ninetails'),
 workspace = system.file('extdata',
                         'test_data',
                         'basecalled_fast5',
                         package = 'ninetails'),
 basecall_group = 'Basecall_1D_000',
 moves = FALSE,
 rescale = TRUE)

print(plot)

} # }
```
