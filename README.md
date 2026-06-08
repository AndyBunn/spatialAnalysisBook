# Spatial Analysis for Environmental Science

Course notes for ESCI 505 at Western Washington University. Started as weekly handouts for a ten-week graduate course, grew into something closer to a short book.

**Live site:** [https://andybunn.github.io/spatialNotes/](https://andybunn.github.io/spatialNotes/)

## What this covers

The notes work through spatial analysis from the ground up, organized into four parts:

- **Foundations** - spatial data structures in R, coordinate systems, methods and generics
- **Point Patterns & Autocorrelation** - point process analysis, Moran's I, LISA
- **Interpolation & Geostatistics** - IDW, thin plate splines, kriging, regression kriging
- **Regression** - GLS, spatial regression with SAR models

The emphasis is on building intuition and getting things done in R. Not on mathematical derivation. We work with point observations and continuous raster surfaces. We simulate known patterns before applying methods to real data, because the best way to trust a tool is to watch it find something you planted.

## Who it's for

Masters students in environmental sciences who need to handle spatial structure in their data. That said, if you're an ecologist or field scientist who wants a practical, R-based entry point to spatial analysis, these notes might work for you too.

## License

This material is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). You're free to use, adapt, and redistribute it for any purpose, including commercial, as long as you give appropriate credit.

## Built with

[Quarto](https://quarto.org), rendered as a book. To build locally you'll need R, Quarto, and the packages loaded in the chapters.
