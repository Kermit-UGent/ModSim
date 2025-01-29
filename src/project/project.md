---
title: "Project description"
tags: ["project"]
order: 1
layout: "md.jlmd"
---

<style>
main a img {
    width: 5rem;
    margin: 1rem;
}
</style>

# Project assignment

**Modelling and Simulation of Biological Systems**

## Assignment

This project aims to apply the principles of the course Modelling and Simulation to a *small*, self-contained example related to bioscience engineering. To this end, you can draw inspiration from our abstracts. You must work on your project in the Pluto notebook environment. 

Your project must be of a maximum length of **eight pages when printed**. Your project contains the following:

- [ ] an abstract with context, why it is relevant, a summary of what you have done and a conclusion
- [ ] a dynamical model (based on ODEs or similar models such as a Markov Chain)
- [ ] a clear outline of the variables and parameters
- [ ] it has some aspect of either process optimization (so you use optimization to improve a parameter, input or decision) or calibration (you tweak a parameter based on sampling or optimization)
- [ ] you use some form of uncertainty assessment, sensitivity analysis or some other stochastic component

These aspects can be extended or limited, as you choose. For example, you use your model for some process optimization via an optimizer but you can also tune a parameter by hand.

You must submit the project as a PDF, HTML, and Julia (.jl) file through UFORA. Between submitting and the final exam, you will be asked to evaluate three other projects in a (limited) manner.

Each project should also end with a small attribution who did what parts according to the [CRediT](https://en.wikipedia.org/wiki/Contributor_Roles_Taxonomy) (Contribution Roles Taxonomy) classification. For example:

> MS: Conceptualization, Methodology, Writing – Review & Editing; BP: Software, visualization, Formal Analysis; DG: Writing – Original Draft Preparation.


## Rubric for grading

| Category | 0-1 Points (Unsatisfactory) | 2 Points (Developing) | 3 Points (Satisfactory) | 4 Points (Good) | 5 Points (Excellent) |
|---|---|---|---|---|---|
| **Clarity and Organization of Notebook** | Notebook is very difficult to follow. Code and explanations are disorganized or missing.  | Notebook has some organization, but it's challenging to understand the logic and purpose. | Notebook follows a generally clear structure with basic explanations of code and results. | Notebook is well-organized, with clear sections and explanations that guide the reader's understanding. | Notebook is exceptionally well-structured, with detailed comments and explanations that make it effortless to follow the project's logic. |
| **Quality of Mathematical Model** | Model is irrelevant to the chosen biological phenomenon or has major conceptual flaws.  | Model shows some relevance to the problem but has significant simplifications or inaccuracies. | Model accurately captures the essential aspects of the biological phenomenon. | Model demonstrates good understanding of the system and includes relevant details and assumptions. | Model is sophisticated and incorporates nuanced or insightful elements that reflect a deep understanding of the biological system. | 
| **Quality of Analysis**  | Analysis is absent or uses incorrect techniques. Results are not presented or interpreted. | Analysis is attempted but flawed (errors, inappropriate methods). Results are presented with limited interpretation. | Analysis uses appropriate techniques and produces mostly correct results with basic interpretation. | Analysis employs suitable techniques leading to correct and meaningful results. Interpretation offers some insights. | Analysis utilizes a range of techniques providing a comprehensive understanding of the model behavior. Interpretation offers significant insights and implications. |

## Tips and advice

- Keep it as simple as possible. A project exploring a small model with two or three variables, outlining a well-known concept from your courses can lead to excellent projects or marks.
- You can use ChatGPT, Gemini or any other AI tool to help you at any point of the project, from brainstorming the initial idea to writing and proofreading text to helping with the code.
- Make sure you explain your reasoning well. Let others read it if it is clear how to do it.
- Make uninteresting or hard-o-understand pieces of code invisible if they don't help the reader.
- You can use `DifferentialEquations.jl` or `ModelingToolkit.jl` if you want to build models that are beyond the scope of `Catalyst.jl`. For example, if you want to model a thermal process. However, this you personal choice, projects that use extra software will not be marked higher. It is perfectly possible to work out an exellent project using only the examples from the practical notes and theory course.
- We expect figures to be tidy (titles, labels and axis). See the example projects for how to do this.
- Use clear variable names and use comments (`#`) to annotate parts that are not clear.
- The programming can and should be minimal.
- Teaching staff will be available for feedback after the lectures and labs or at designated time slots.
- Keep it as simple as possible. 