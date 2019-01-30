# Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)

function info()
    println("Package name: ClinicalTrialUtilities")
    println("Clinical Trial Power and Sample Size Utilities")
    println("Copyright © 2019 Vladimir Arnautov aka PharmCat (mail@pharmcat.net)")
    println("For references, citations and copyrights type <ClinicalTrialUtilities.citation()>")
    println("Licence <ClinicalTrialUtilities.licence()>")
end

function citation()
    println("References:")
    println("")
    println("PowerTOST Authors©R :Detlew Labes (DetlewLabes@gmx.de), Helmut Schuetz (helmut.schuetz@bebac.at), Benjamin Lang (lang.be@icloud.com)")
    println("")
    println("Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series.")
    println("")
    println("Connor R. J. 1987. Sample size for testing differences in proportions for the paired-sample design. Biometrics 43(1):207-211. page 209.")
    println("")
    println("Owen, D B (1965) \"A Special Case of a Bivariate Non-central t-Distribution\" Biometrika Vol. 52, pp.437-446.")
    println("")
    println("D.B. Owen \"Tables for computing bivariate normal Probabilities\" The Annals of Mathematical Statistics, Vol. 27 (4) Dec. 1956, pp. 1075-1090")
    println("")
    println("FORTRAN code by J. Burkhardt, license GNU LGPL")
    println("")
    println("Matlab code  by J. Burkhardt license GNU LGPL")
    println("")
end

function licence()
    println("Licence:")
    println("GNU Affero General Public License v3.0")
    println("see <https://www.gnu.org/licenses/>")
end
