# Spatial Analysis for Environmental Data

*An R-Based Introduction*

**Read it here: [spatial.andybunn.org](https://spatial.andybunn.org/)**

This started as weekly handouts for ESCI 505, a graduate course in spatial analysis at Western Washington University, and grew into a short book. It's free to read online and openly licensed.

## What it covers

The book works through spatial analysis from the ground up, in four parts:

- **Foundations** - spatial data structures in R, coordinate systems, methods and generics
- **Point Patterns & Autocorrelation** - point process analysis, Moran's I, LISA
- **Interpolation & Geostatistics** - IDW, thin plate splines, kriging, regression kriging
- **Regression** - GLS, spatial regression with SAR models

The emphasis is on building intuition and getting things done in R, not on mathematical derivation. We stick to point observations and continuous raster surfaces, and we simulate known patterns before turning a method loose on real data, because the best way to trust a tool is to watch it find something you planted.

## Who it's for

Orignally? Masters students in environmental science who need to handle spatial structure in their data. But if you're an ecologist or field scientist who wants a practical, R-based entry point to spatial analysis, it should work for you too. It assumes an introductory linear-modeling course and not math beyond calculus (and doesn't even need that really).

## How it's built

The book is written in [Quarto](https://quarto.org) and rendered as a book site. Every push to `main` triggers a GitHub Actions workflow that renders the book and deploys it to [spatial.andybunn.org](https://spatial.andybunn.org/), so the live site is always current.

To build it locally you'll need R, Quarto, and the packages listed in `DESCRIPTION`. From the project root:

```sh
quarto render
```

Found an error or have a suggestion? Open an issue, or use the "Edit this page" link on any chapter. I'd love to hear from you.

## License

The book is dual-licensed (see [LICENSE.md](LICENSE.md)):

- **Prose and figures** under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/). Share and adapt for noncommercial use, keep derivatives open under the same terms, and give credit.
- **Code** under the [MIT License](https://opensource.org/licenses/MIT). Use it however you like.
