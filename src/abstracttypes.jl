abstract type AbstractData end
abstract type AbstractDataSet end

abstract type AbstractSubject <: AbstractData end


abstract type AbstractTask end
abstract type AbstractParameter end

abstract type AbstractProportion  <:  AbstractParameter end
abstract type AbstractMean  <:  AbstractParameter end

abstract type AbstractSimpleProportion <:  AbstractProportion end

abstract type AbstractCompositeProportion  <:  AbstractProportion end
abstract type AbstractCompositeMean{T}  <:  AbstractMean end

abstract type AbstractObjective end
abstract type AbstractSampleSize <: AbstractObjective end
abstract type AbstractPower <: AbstractObjective end

abstract type AbstractHypothesis end
abstract type AbstractEquivalenceHypothesis <: AbstractHypothesis end

abstract type AbstractDesign end
