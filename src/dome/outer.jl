
"""
The outer dome is the intersection of the triangles formed by the upper left and right legs of the corridors along with the x axis.
We compute it by finding the lowest intersection of a left leg with a right leg, then finding the fastest way to get to y=0 from there.
We compute the intersection points ordered by x coordinate, starting from the leftmost point.
"""
function outerDome( outerLRays::SubArray{LRay, 1}, outerTs::Array{Rational, 4} ) :: Tuple{Matrix{Rational}, Array{Integer,1}, Array{Integer,1}}
    xPeak, yPeak, leftIndex, rightIndex = findOuterDomePeak( outerLRays, outerTs )
 
    @assert yPeak != Inf # a peak should always be found
    if yPeak <=0  return Matrix{Rational}(undef, 0, 0), [], []  end # lowest intersection is not above x axis

    # compute intersection times for x-axis
    xAxisIntersectionTs = getIntersectionTsWithXAxis( outerLRays )

    # compute left and right indices of outer dome
    outerDomeIndicesL = outerDomeIndicesFromPeakLeft(  outerLRays, outerTs, leftIndex, rightIndex, @view xAxisIntersectionTs[:,2,LEFT] )
    outerDomeIndicesR = outerDomeIndicesFromPeakRight( outerLRays, outerTs, leftIndex, rightIndex, @view xAxisIntersectionTs[:,2,RIGHT] )

    # compute number of points on the outer dome
    nL, nR = length( neighbors( outerDomeIndicesL ) ), length( neighbors( outerDomeIndicesR ) )
    n = 1 + nL + 1 + nR + 1 # +1 for each intersection with x axis, +1 for the peak 

    # matrix, where every column is a point on the outer dome, ordered by x coordinate
    outerDome = zeros( Rational, 2, n )

    # compute the intersection points with x-axis (y=0)
    firstT, lastT = xAxisIntersectionTs[outerDomeIndicesL[1],1,LEFT], xAxisIntersectionTs[outerDomeIndicesR[end],1,RIGHT] # [index, time wrt x-axis, side]
    outerDome[1,1]   = firstT # we conveniently chose the x axis intersection t to be the x value of the intersection point
    outerDome[1,end] = lastT  # that's why we don't need to step the ray to get the x value

    # compute remaining points on the outer dome
    if nL>0
        outerDome[:,2:nL+1]     = reduce( hcat, [ stepRayByT( leftRay(  outerLRays[i] ), outerTs[i,j,LEFT,LEFT] )   for (i,j) in neighbors( outerDomeIndicesL ) ] )
    end
    outerDome[:,nL+2]           = [ xPeak, yPeak ]
    if nR>0
        outerDome[:,nL+3:end-1] = reduce( hcat, [ stepRayByT( rightRay( outerLRays[i] ), outerTs[i,j,RIGHT,RIGHT] ) for (i,j) in neighbors( outerDomeIndicesR ) ] )
    end

    return outerDome, outerDomeIndicesL, outerDomeIndicesR
end

# function that computes the truncated outer dome, i.e. the outer dome with everything above y=1 removed
function outerDomeTruncated( outerLRays::SubArray{LRay, 1}, outerTs::Array{Rational, 4} ) :: Tuple{Matrix{Rational}, Array{Integer,1}, Array{Integer,1}}
    oDome, outerDomeIndicesL, outerDomeIndicesR = outerDome( outerLRays, outerTs )
    if isempty( oDome ) return oDome, outerDomeIndicesL, outerDomeIndicesR end
    return truncate( oDome ), outerDomeIndicesL, outerDomeIndicesR
end

# finds the lowest intersection of a left ray with a right ray
# assumes no two rays are parallel
function findOuterDomePeak( outerLRays::SubArray{LRay, 1}, outerTs :: Array{Rational, 4} ) :: Tuple{Rational, Rational, Integer, Integer}
    k = length( outerLRays )

    # initialize minimum y coordinate, indices of LRays that intersect at that point
    # leftIndex is the index of the LRays that intersect at the point with its left leg
    # rightIndex is the index of the LRays that intersect at the point with its right leg
    xMin, yMin, leftIndex, rightIndex = Inf, Inf, 0, 0

    # first find minimal intersection of left rays with right rays on same LRays (corner points)
    for i in 1:k
        local lRay = outerLRays[i]
        if lRay.origin[2] < yMin
            xMin, yMin = lRay.origin
            leftIndex = rightIndex = i
            if lRay.origin[2] <= 0 # not above x axis
                return xMin, yMin, leftIndex, rightIndex
            end
        end
    end

    # then check the other pairs of rays
    for i in 1:k
        for j in i+1:k
            local isOnIthLegLeft  = outerTs[i,j,LEFT ,RIGHT] <= 0
            local isOnIthLegRight = outerTs[i,j,RIGHT,LEFT] >= 0
            local isOnJthLegLeft  = outerTs[j,i,LEFT ,RIGHT] <= 0
            local isOnJthLegRight = outerTs[j,i,RIGHT,LEFT] >= 0
            if isOnIthLegLeft && isOnJthLegRight
                xLR, yLR = stepRayByT( rightRay( outerLRays[j] ), outerTs[j,i,RIGHT,LEFT] )
                if yLR < yMin
                    xMin, yMin = xLR, yLR
                    leftIndex, rightIndex = i, j
                elseif yLR == yMin # if y is equal to yMin, then we need to check if other legs are steeper (and therefore make the dome smaller)
                    leftIndex  = tiebreakers[LEFT, Integer( FAST ) + 1]( [i,leftIndex],  outerLRays )
                    rightIndex = tiebreakers[RIGHT,Integer( FAST ) + 1]( [i,rightIndex], outerLRays )
                end
            end
            if isOnIthLegRight && isOnJthLegLeft
                xRL, yRL = stepRayByT( rightRay( outerLRays[i] ), outerTs[i,j,RIGHT,LEFT] )
                if yRL < yMin
                    xMin, yMin = xRL, yRL
                    leftIndex, rightIndex = j, i
                elseif yRL == yMin # if y is equal to yMin, then we need to check if other legs are steeper (and therefore make the dome smaller)
                    leftIndex  = tiebreakers[LEFT, Integer( FAST ) + 1]( [j,leftIndex],  outerLRays )
                    rightIndex = tiebreakers[RIGHT,Integer( FAST ) + 1]( [i,rightIndex], outerLRays )
                end
            end
        end
    end
    return xMin, yMin, leftIndex, rightIndex
end

"""
Compute outer dome indices from a peak on the right side.
"""
outerDomeIndicesFromPeakRight( outerLRays::SubArray{LRay,1}, outerTs::Array{Rational,4}, leftIndex::Integer, rightIndex::Integer, xAxisIntersectionTs::SubArray{Rational,1} ) :: Array{Integer,1} =
    descend( outerLRays, outerTs, leftIndex, rightIndex, xAxisIntersectionTs, RIGHT, FAST )[1:end-1] # remove last index to get rid of 0

"""
Compute outer dome indices from a peak on the left side.
"""
outerDomeIndicesFromPeakLeft( outerLRays::SubArray{LRay,1}, outerTs::Array{Rational,4}, leftIndex::Integer, rightIndex::Integer, xAxisIntersectionTs::SubArray{Rational,1}  ) :: Array{Integer,1} =
    # reverse indices to get them in order of increasing x coordinate
    reverse!( descend( outerLRays, outerTs, rightIndex, leftIndex, xAxisIntersectionTs, LEFT, FAST )[1:end-1] ) # remove last index to get rid of 0