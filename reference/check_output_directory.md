# Check and handle existing output directory for ninetails analysis

This function checks if the specified output directory already exists
and contains files that might be overwritten by the ninetails analysis.
If the directory exists and is not empty, it prompts the user for action
and logs the decision.

## Usage

``` r
check_output_directory(save_dir, log_message, input_fn = readline)
```

## Arguments

- save_dir:

  Character string. Full path to the output directory where ninetails
  results will be saved.

- log_message:

  Function for logging messages to both console and log file. Should
  accept parameters: message, type, section.

- input_fn:

  Function for reading user input. Defaults to `readline`. Can be
  overridden for testing purposes.

## Value

Logical. Returns TRUE if the analysis should proceed, FALSE if the user
chose to abort the analysis.

## Details

This function is not intended to be used outside the pipeline wrapper.

## User Interaction

When an existing non-empty directory is detected, the function will:

- Display the directory path and file count

- Prompt the user to choose: abort analysis or overwrite existing files

- Wait for user input (a/A for abort, o/O for overwrite)

- Log the user's decision and proceed accordingly

## Directory States

The function handles several directory states:

- **Non-existent**: Creates the directory and logs creation

- **Empty**: Uses existing directory and logs confirmation

- **Non-empty**: Prompts user for overwrite decision

## Implementation Notes

- Uses [`readline()`](https://rdrr.io/r/base/readline.html) for
  interactive user input

- Logs all decisions and actions for audit trail

- Handles edge cases like permission errors

- Validates user input with retry mechanism

## Examples

``` r
if (FALSE) { # \dontrun{
# Example with logging function
log_func <- function(msg, type = "INFO", section = NULL) {
  cat(sprintf("[%s] %s\n", type, msg))
}

# Check directory and proceed if allowed
should_proceed <- check_output_directory("/path/to/output", log_func)
if (should_proceed) {
  # Continue with analysis
} else {
  # User chose to abort
}
} # }
```
