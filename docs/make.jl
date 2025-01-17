using Documenter, ClinicalTrialUtilities, DataFrames

makedocs(
    modules = [ClinicalTrialUtilities],
    sitename = "ClinicalTrialUtilities",
    authors = "Vladimir Arnautov",
    linkcheck = false,
    doctest = false,
    pages = [
        "Home" => "index.md",
        "Sample size" => "samplesize.md",
        "Power" => "power.md",
        "Randomization" => "random.md",
        "Utilities" => "utils.md",
        "Examples" => "examples.md",
        "References" => "ref.md",
    ],
)

deploydocs(repo = "github.com/PharmCat/ClinicalTrialUtilities.jl.git", push_preview = true)
