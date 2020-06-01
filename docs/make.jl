using Documenter, ClinicalTrialUtilities, DataFrames, CSV, Plots

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
        "Confidence intervals" => "ci.md",
        "Descriptive statistics" => "ds.md",
        "NCA" => ["NCA" => "nca.md",
        "Pharmakokinetics" => "pk.md",
        "Pharmacodynamics" => "pd.md",
        "PK Plots" => "pkplot.md"],
        "Randomization" => "random.md",
        "Simulations" => "sim.md",
        "Utilities" => "utils.md",
        "Examples" => "examples.md",
        "Export" => "export.md",
        "Validation" => "validation.md",
        "References" => "ref.md",
    ],
)

deploydocs(repo = "github.com/PharmCat/ClinicalTrialUtilities.jl.git")
