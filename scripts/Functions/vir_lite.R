
# From https://stackoverflow.com/questions/49166411/set-transparency-saturation-of-palette-in-ggplot

vir_lite <- function(cols, ds=0.5, dv=0.8) {
  cols = rgb2hsv(col2rgb(cols))
  cols["v", ] = cols["v", ] + dv*(1 - cols["v", ])
  cols["s", ] = ds*cols["s", ]
  apply(cols, 2, function(x) hsv(x[1], x[2], x[3]))
}
