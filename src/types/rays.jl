
struct Ray
    origin    :: Point2d
    direction :: Point2d
end

struct LRay
    origin :: Point2d
    left   :: Point2d
    right  :: Point2d
    LRay( origin::Point2d, left::Point2d ) = new( origin, left, rotate90CW( left ) )
end

rotate( lRay::LRay, triple::PythagoreanTriple ) :: LRay =
    LRay( rotate( lRay.origin, triple ), rotate( lRay.left, triple ) )

function lRaysFromCornersAndTriple( innerCorner, outerCorner, triple ) :: Array{LRay, 1}
    # rotate inner and outer corner by α
    innerOrigin = rotate( innerCorner, triple )
    outerOrigin = rotate( outerCorner, triple )

    # compute the direction of the left leg
    leftLegDirection = rotate( HORIZONTAL, triple )

    # create the inner and outer LRay
    inner = LRay( innerOrigin, leftLegDirection )
    outer = LRay( outerOrigin, leftLegDirection )

    return [ inner, outer ]
end

function getIntersectionTs( r1::Ray, r2::Ray ) :: SVector{2,Rational}
    # Define matrix A and vector b based on the equations of the rays
    A = hcat( r1.direction, -r2.direction )
    b = safely( -, r2.origin, r1.origin )

    # Solve the system of linear equations using \ operator
    t = safely( \, A, b )
    return t
end

intersectionPointFromRays( r1::Ray, r2::Ray ) :: Point2d = stepRayByT( r1, getIntersectionTs( r1, r2 )[1] )

function getIntersectionTs( lray1::LRay, lray2::LRay ) :: SizedArray{Tuple{2,2,2,2},Rational}
    ts = SizedArray{Tuple{2,2,2,2},Rational}( zeros( Rational, 2, 2, 2, 2 ) )

    r1L, r1R = leftRay( lray1 ), rightRay( lray1 )
    r2L, r2R = leftRay( lray2 ), rightRay( lray2 )

    # ray 1 left with ray 2 left
    ts[1, 2, 1, 1], ts[2, 1, 1, 1] = getIntersectionTs( r1L, r2L )

    # ray 1 left with ray 2 right
    ts[1, 2, 1, 2], ts[2, 1, 1, 2] = getIntersectionTs( r1L, r2R )

    # ray 1 right with ray 2 left
    ts[1, 2, 2, 1], ts[2, 1, 2, 1] = getIntersectionTs( r1R, r2L )

    # ray 1 right with ray 2 right
    ts[1, 2, 2, 2], ts[2, 1, 2, 2] = getIntersectionTs( r1R, r2R )

    return ts
end


function getIntersectionTs( rs::SubArray{LRay,1} ) :: Array{Rational,4}
    n = length( rs )
    # Create 4D array to store the results
    # n LRays x n LRays x 2 directions x 2 directions, t is always in the first ray's direction
    # we initialize the array with zeros, since we know that the intersection of a ray with itself or its opposite is always 0
    ts = zeros( Rational, n, n, 2, 2 )

    # Loop over all pairs of rays
    for i in 1:n
        for j in i+1:n
            ts_ij = getIntersectionTs( rs[i], rs[j] )
            ts[i, j, :, :] = ts_ij[1, 2, :, :]
            ts[j, i, 1, :] = ts_ij[2, 1, :, 1]
            ts[j, i, 2, :] = ts_ij[2, 1, :, 2]
        end
    end

    return ts
end

"""
    getIntersectionTsWithXAxis( r::Ray ) :: SVector{2,Rational}

    Compute the intersection times of a Ray with the x-axis.
    Returns an array with two elements:
        1st element: intersection time wrt to x-axis
        2nd element: intersection time wrt to ray's direction
"""
function getIntersectionTsWithXAxis( r::Ray ) :: SVector{2,Rational}
    try
        # equivalent to `hcat( HORIZONTAL, -r.direction ) \ r.origin` but faster
        yDivY = r.origin[2]//r.direction[2]
        return SA[ r.origin[1] - r.direction[1]*yDivY, -yDivY ]
    catch e
        if isa( e, OverflowError )
            yDivY = big(r.origin[2])//big(r.direction[2])
            return SA[ big(r.origin[1]) - big(r.direction[1])*yDivY, -yDivY ]
        else
            rethrow( e )
        end
    end
end

"""
    getIntersectionTsWithXAxis( lRay::LRay ) :: SMatrix{2,2,Rational}

    Compute the intersection times of a Ray with the x-axis.
    Returns an array with two elements:
        1st element: intersection time wrt to x-axis
        2nd element: intersection time wrt to ray's direction
"""
function getIntersectionTsWithXAxis( lRay::LRay ) :: SMatrix{2,2,Rational}
    solLeft  = getIntersectionTsWithXAxis( leftRay( lRay ) )
    solRight = getIntersectionTsWithXAxis( rightRay( lRay ) )
    return hcat( solLeft, solRight )
end

"""
    getIntersectionTsWithXAxis( lRays::SubArray{LRay,1} ) :: Array{Rational,3}

    Compute the intersection times of a set of LRays with the x-axis.
    The dimensions of the returned array have the following meaning:
        1st dimension: index of the LRay
        2nd dimension: intersection time wrt to (1: x-axis, 2: ray's direction)
        3rd dimension: side of the corridor (1: left, 2: right)
"""
getIntersectionTsWithXAxis( lRays::SubArray{LRay,1} ) :: Array{Rational,3} =
 permutedims( cat( dims=3, getIntersectionTsWithXAxis.( lRays )... ), (3,1,2) )

stepRayByT( r::Ray, t::Rational ) :: Point2d = try 
    r.origin + t * r.direction
    catch e
        if isa( e, OverflowError )
            big(r.origin) + big(t) * big(r.direction)
        else
            rethrow( e )
        end
    end
stepLRayByT( lRay::LRay, side::Integer, t::Rational ) :: Point2d = stepRayByT( side == LEFT ? leftRay( lRay ) : rightRay( lRay ), t )

rightRay( lRay::LRay ) :: Ray = Ray( lRay.origin, lRay.right )
leftRay( lRay::LRay )  :: Ray = Ray( lRay.origin, lRay.left )

 # slope computing
slope( p::Point2d ) :: Rational = p[2]//p[1]
slope( r::Ray ) :: Rational = slope( r.direction )
slope( lRay::LRay, side::Integer ) :: Rational = slope( side==1 ? lRay.left : lRay.right ) 

function minkowskiLRaysFromBoxAndTriples( box::SubBox, triples::TripleArray ) :: Array{LRay,2}
    k = length( triples )
    lRays = Array{LRay,2}(undef, k, 2)
    for i in 1:k
        triple = triples[i]

        xlo, xhi = boundsFromSubInterval( box[2i-1] )
        ylo, yhi = boundsFromSubInterval( box[2i] )
        innerCorner = SA[xlo, ylo]
        outerCorner = SA[xhi + 1, yhi + 1]

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