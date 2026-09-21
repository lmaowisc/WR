# Plot the standardized score processes

Plot the standardized score processes.

## Usage

``` r
# S3 method for class 'pwreg.score'
plot(
  x,
  k,
  xlab = "Time",
  ylab = "Standardized score",
  lty = 1,
  frame.plot = TRUE,
  add = FALSE,
  ylim = c(-3, 3),
  xlim = NULL,
  lwd = 1,
  ...
)
```

## Arguments

- x:

  an object of class `pwreg.score`.

- k:

  A positive integer indicating the order of covariate to be plotted.
  For example, `k=3` requests the standardized score process for the
  third covariate in the covariate matrix `Z`.

- xlab:

  a title for the x axis.

- ylab:

  a title for the y axis.

- lty:

  the line type. Default is 1.

- frame.plot:

  a logical variable indicating if a frame should be drawn in the 1D
  case.

- add:

  a logical variable indicating whether add to current plot?

- ylim:

  a vector indicating the range of y-axis. Default is (-3,3).

- xlim:

  a vector indicating the range of x-axis. Default is NULL.

- lwd:

  the line width, a positive number. Default is 1.

- ...:

  further arguments passed to or from other methods.

## Value

A plot of the standardized score process for object `pwreg.score`.

## See also

[`score.proc`](https://lmaowisc.github.io/WR/reference/score.proc.md)

## Examples

``` r
# see the example for score.proc
```
