module Sofa

export Point2d, Interval, SubInterval
export Ray, LRay, Corridor, PythagoreanTriple
export getIntersectionTs, getIntersectionTsWithXAxis
export rotate, rotate90CW
export stepRayByT, outerDome, innerDome, findOuterDomePeak
export lRaysFromCorridor, lRaysFromCorridors
export truncate, innerDomeTruncated, outerDomeTruncated
export getCorridorIntersectionAreas
export initialBoxFromTriples, subBoxFromTriples
export midpoint, refineBounds
export computeBounds
export genTriples, intermediateTriples, isPythagoreanTriple
export minkowskiAreaFromBoxAndTriples, minkowskiAreaFromBoxesAndTriples
export minDome, maxDome
export sampleCorridorsFromTriples, sampleTotalCorridorIntersectionArea, isInsideCorridors

# export plotPoints!, plotRay!, plotLRay!, plotLRays! # not part of package, add manually if needed

using LinearAlgebra
using StaticArrays


include("utils.jl")

include("types/point2d.jl")
include("types/pythagoreanTriple.jl")
include("types/interval.jl")
include("types/box.jl")
include("types/rays.jl")
include("types/corridor.jl")

include("dome/helpers.jl")
include("dome/outer.jl")
include("dome/inner.jl")

include("bounds.jl")
include("compute.jl")
include("stateTools.jl")

# include("plotting/plotting.jl") # not part of package, add manually if needed

end