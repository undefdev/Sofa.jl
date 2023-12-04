
# Corridor is a struct that represents a corridor in the plane. It is defined by the inner corner and the angle of rotation,
# represented by a Pythagorean triple (a, b, c) where a² + b² = c².
# The corridor is the area enclosed by the two LRays that are created from the inner corner and rotated by the angle of rotation.
# The outer LRay is the inner LRay translated by (1,1) before rotation.
struct Corridor
    innerCorner::Point2d # point in ℝ²
    rotation::PythagoreanTriple # angle of rotation ∈ (0, π/2), represented by a Pythagorean triple
end

function isInsideCorridor( p::Point2d, corridor::Corridor ) :: Bool
    # rotate p by -α
    rot = (-corridor.rotation[1], corridor.rotation[2], corridor.rotation[3])
    pRotated = rotate( p, rot )
    # translate pRotated by -innerCorner
    pTranslated = pRotated - corridor.innerCorner
    # check if pRotated is inside standard corridor
    return 0 <= maximum( pTranslated ) <= 1
end

isInsideCorridors( p::Point2d, corridors::Array{Corridor,1} ) :: Bool = all( isInsideCorridor( p, corridor ) for corridor in corridors )

function lRaysFromCorridor( corridor::Corridor ) :: Array{LRay, 1}
    return lRaysFromCornersAndTriple( corridor.innerCorner, corridor.innerCorner + [1, 1], corridor.rotation )
end

# returns a k×2 matrix of LRays, where the first column contains the inner LRays and the second column the outer LRays
function lRaysFromCorridors( corridors::Array{Corridor,1} ) :: Matrix{LRay}
    return permutedims( hcat( lRaysFromCorridor.( corridors )... ), (2,1) )
end

function sampleCorridorsFromBoxAndTriples( box::Box, triples::TripleArray ) :: Array{Corridor,1}
    local k = length( triples )
    local corridors = Array{Corridor,1}(undef, k)
    for i in 1:k
        local tx, ty = rand(1:100, 2).//100
        local x, y = lerp( box[2i-1], tx ), lerp( box[2i], ty )
        corridors[i] = Corridor( SA[ x, y ], triples[i] )
    end
    return corridors
end

function sampleCorridorsFromTriples( triples::TripleArray ) :: Array{Corridor,1}
    return sampleCorridorsFromBoxAndTriples( initialBoxFromTriples( triples ), triples )
end

function getCorridorIntersectionAreas( corridors::Array{Corridor,1} ) :: Array{Rational,1}
    # get outer and inner LRays from corridors
    lRays = lRaysFromCorridors( corridors )
    return areasBetweenLRays( lRays )
end

function areasBetweenLRays( lRays::Array{LRay,2} ) :: Array{Rational{BigInt},1}
    innerLRays = @view lRays[:,1]
    outerLRays = @view lRays[:,2]

    # get intersection times of inner and outer LRays
    innerTs, outerTs = getIntersectionTs( innerLRays ), getIntersectionTs( outerLRays )

    # compute domes
    iDomeComponentsTruncated, _ = innerDomeTruncated( innerLRays, innerTs )
    oDomeTruncated, _, _        = outerDomeTruncated( outerLRays, outerTs )

    if isempty( oDomeTruncated ) return [] end # no intersection

    # compute piecewise linear function representing the difference of the outer dome and the inner domes (where each dome is seen as piecewise linear function)
    heightDifferences = subtractPiecewise( oDomeTruncated, iDomeComponentsTruncated )

    # we only care about the area between the domes when the outer dome is above the inner dome, so we set all negative y-values to 0
    heightDifferences[2, :] = max.(heightDifferences[2, :], 0)

    # we want the area of the connected components of the intersection, so we want to split the total integrand by its support intervals

    # find all indices of zero y-values of the total integrand
    zeroIndices = findall( heightDifferences[2,:] .== 0 )

    # get indices to split the total integrand by its support intervals
    componentRanges = [ iZero:iNextZero for (iZero, iNextZero) in neighbors( zeroIndices ) if iNextZero - iZero > 1 ]
    
    # integrate views of the total integrand and return
    return [ compute_area( @view heightDifferences[:,componentRange] ) for componentRange in componentRanges ]
end

function maxAreaBetweenLRays( lRays::Array{LRay,2} ) :: Rational{BigInt}
    integrals = areasBetweenLRays( lRays )
    return isempty( integrals ) ? 0//1 : maximum( integrals )
end

function centerLRaysFromBoxAndTriples( box::SubBox, triples::TripleArray ) :: Array{LRay,2}
    k = length( triples )
    lRays = Array{LRay,2}(undef, k, 2)
    midpoint = mid.(box)
    for i in 1:k
        triple = triples[i]

        x, y = midpoint[2i-1], midpoint[2i]
        innerCorner = SA[x, y]
        outerCorner = innerCorner .+ 1

        # rotate inner and outer corner by α
        innerOrigin = rotate( innerCorner, triple )
        outerOrigin = rotate( outerCorner, triple )

        # compute the direction of the left leg
        leftLegDirection = rotate( HORIZONTAL, triple )

        # create the inner and outer LRay
        lRays[i,1] = LRay( innerOrigin, leftLegDirection )
        lRays[i,2] = LRay( outerOrigin, leftLegDirection )
    end
    return lRays
end

function centerAreaFromBoxAndTriples( box::SubBox, triples::TripleArray ) :: Rational{BigInt}
    lRays = centerLRaysFromBoxAndTriples( box, triples )
    return maxAreaBetweenLRays( lRays )
end

function minkowskiAreaFromBoxAndTriples( box::SubBox, triples::TripleArray ) :: Rational{BigInt}
    lRays = minkowskiLRaysFromBoxAndTriples( box, triples )
    return maxAreaBetweenLRays( lRays )
end

using Distributions
# function to compute the total area of the intersection of a set of corridors in ℝ² intersected with H:=ℝ×[0,1]
function sampleTotalCorridorIntersectionArea( corridors::Array{Corridor,1}, nSamples::Integer, boundXLeft::Real=-6, boundXRight::Real=6 ) :: Real
    dx, dy = Uniform( boundXLeft, boundXRight ), Uniform( 0, 1 )
    # sample points in [boundXLeft, boundXRight]×[0,1] and check if they are inside the corridors and H
    points = [ SA[ rand( dx ), rand( dy ) ] for i in 1:nSamples ]
    return sum( isInsideCorridors( p, corridors ) for p in points ) / nSamples * ( boundXRight - boundXLeft )
end
