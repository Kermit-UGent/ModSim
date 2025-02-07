---
title: "Introduction"
tags: ["exercises"]
order: 1
layout: "md.jlmd"
---

<style>
main a img {
    width: 5rem;
    margin: 1rem;
}
</style>

# Exercises description

Here are renders of the exercises, see all pages on the left.

To download all exercises: see Ufora.

We have annotated the exercises with either a number or an extra prefix.

> 1. XYZ

Will be covered in exercise lession 1.

> EXTRA. XYZ

Additional exercises that will not be covered in the guided exercises.

# Notes on the dependencies

If you insist on downloading the exercises from this website, note that because we are rendering the notebooks here, we make use of our a specific environment. You will need to update this on your system. Look out for the cell with:
```julia
using Pkg
Pkg.activate("../../pluto-deployment-environment")
```
Change this to the your current folder so that the Project and Manifest files are generated there:
```julia
using Pkg
Pkg.activate(".") 
```
