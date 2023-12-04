
# Constants to improve code readability
const LEFT, RIGHT = 1, 2
@enum DescendMode SLOW FAST

"""
`bestSlope( candidateIndices::Array{<:Integer,1}, lRays::SubArray{LRay,1}, side::Integer, mode::DescendMode ) :: Integer`

Finds the "best" index from a list of candidates based on the slope of the corresponding ray.
What "best" means depends on the side and mode.
The left leg of a ray has a positive slope, the right leg has a negative slope.
If mode is SLOW, we want to descend as slowly as possible, so for a left leg we want small slopes and for a right leg we want our slopes to have small absolute value.
If mode is FAST, we want to descend as quickly as possible, so for a left leg we want large slopes and for a right leg we want our slopes to have large absolute value.
"""
function bestSlope( candidateIndices::Array{<:Integer,1}, lRays::SubArray{LRay,1}, side::Integer, mode::DescendMode ) :: Integer
    (argoptL, argoptR) = mode == FAST ? (argmax, argmin) : (argmin, argmax)
    argopt = side == LEFT ? argoptL : argoptR

    # find the index with that value that corresponds to the leg with the largest slope and return it
    return argopt( i -> slope( lRays[i], side ), candidateIndices )
end

const tiebreakers = [ (cs, lRays) -> bestSlope( cs, lRays, side, mode ) for side in (LEFT, RIGHT), mode in (SLOW, FAST) ]

"""
    descend( lRays::SubArray{LRay,1}, ts::Array{Rational,4}, firstIndex::Integer, secondIndex::Integer, side::Integer, mode::DescendMode ) :: Array{Integer,1}

    Descend along the legs of the rays in `lRays` starting from `firstIndex` and `secondIndex` until the x axis is reached, or additionally in case of SLOW mode, until a leg on the other side is reached.
    The side of the legs to descend on is given by `side`.
    The mode determines how the descent is performed (FAST or SLOW).
    The function returns the indices of the legs that were descended on.
    The value 0 is appended to the list of indices to indicate that the x axis was reached.
"""
function descend( 
    lRays::SubArray{LRay,1},
    ts::Array{Rational,4}, 
    firstIndex::Integer, 
    secondIndex::Integer, 
    xAxisIntersectionTs::SubArray{Rational,1},
    side::Integer,
    mode::DescendMode, 
) :: Array{Integer,1}
    # Initialize indices with first and second index
    indices = [ firstIndex, secondIndex ]

    # Define side specific functions
    local afterT = side == RIGHT ? ((a,t) -> a>t>=0) : ((a,t) -> a<t<=0)
    local afterOrEqT = side == RIGHT ? ((a,t) -> a>=t>0) : ((a,t) -> a<=t<0)
    local closest = side == RIGHT ? minimum : maximum
    local tiebreaker = tiebreakers[side, Integer( mode ) + 1]

    while true
        # Take the last two indices
        local iLast, i = indices[end-1:end]

        local otherSide = 3 - side
        
        # Determine the side based on the value of iLast, every intersection except the first one is between legs on the same side
        lastSide = length(indices)==2 ? otherSide : side
        
        # Get intersection t in terms of the current index
        t = ts[i,iLast,side,lastSide]
        
        # Find candidates for the next index based on their intersection with the ray
        candidateIndices = afterT.( ts[i,:,side,side], t )

        # If the mode is SLOW, we need to remove candidates that don't really intersect with the ray
        # We do this by checking if the intersection time wrt. the candidate is
        # - <0 for right legs (since time is positive on the right leg)
        # - >0 for left legs (since time is negative on the left leg)
        if mode == SLOW
            candidateIndices = candidateIndices .& ( side == RIGHT ? ts[:,i,side,side] .>= 0 : ts[:,i,side,side] .<= 0 )
        end
        
        # If there aren't any candidates, return indices without the first one        
        if !any( candidateIndices )  break  end # we've reached the x-axis. push 0 and return indices[2:end]

        # Get the next t based on the side (for left t<=0, for right t>=0)
        tNext = closest( ts[i,candidateIndices,side,side] )

        # If the mode is SLOW, we want to check if there is an intersection with a leg on the other side no further than tNext, if so we stop descending
        if mode == SLOW
            # Get candidate indices for the other side
            onOtherLeg = otherSide==RIGHT ? >(0) : <(0)
            otherCandidateIndices = afterT.( ts[i,:,side,3-side], t ) .& onOtherLeg.( ts[:,i,3-side,side] )

            if any( otherCandidateIndices )
                # Get the next t for intersections with other sided legs
                tNextOtherSide = closest( ts[i,otherCandidateIndices,side,3-side] )
          
                # If there is a closer or equally close intersection with a leg on the other side, we stop descending
                if afterOrEqT( tNext, tNextOtherSide )
                    # If the next intersection is below or on the x axis, we stop descending and append 0 to the indices
                    if afterOrEqT( tNextOtherSide, xAxisIntersectionTs[i] )  break  end
                    # We didn't reach the x axis, so we don't append 0 to the indices
                    # Return indices without the first one, since it's the same as the second one
                    return indices[2:end]
                end
            end
        end

        # If the next intersection is below or on the x axis, we stop descending and append 0 to the indices
        if afterOrEqT( tNext, xAxisIntersectionTs[i] )  break  end

        # Find candidate indices equal to the next t
        candidateIndices = ts[i,:,side,side] .== tNext

        # Find index with steepest slope
        iNext = tiebreaker( findall( candidateIndices ), lRays )

        # Add the next index to the list
        push!( indices, iNext )
    end

    # We append 0 to indices to indicate that we've reached the x axis.
    push!( indices, 0 )
    return indices[2:end] # return indices without the first one, since it's the same as the second one
 end

# function to truncate a polygon represented by a matrix of vertices, such that everything above y=1 is removed
function truncate( polygon::Matrix{<:Rational} ) :: Matrix{<:Rational}
    # if there is no vertex above y=1, return the polygon
    first = findfirst( >(1), polygon[2,:] )
    if first === nothing  return polygon  end
    # initialize matrix with vertices before the ith vertex
    truncated :: Matrix{<:Rational} = polygon[:,1:first-1]

    # define the line at y=1
    yOne = Ray( SA[0//1, 1//1], HORIZONTAL )
    # iterate over the vertices of the polygon
    for i in first:size( polygon, 2 )
        # if the vertex is below y=1, add it to the truncated polygon
        if polygon[2, i] <= 1
            truncated = [truncated polygon[:, i]]
        else
            # if the vertex is above y=1 and the previous vertex is below y=1, find the intersection point and add it to the truncated polygon
            if polygon[2, i-1] < 1 # we know i>1, since the first vertex is on the x axis
                ray = Ray( polygon[:, i-1], polygon[:, i] - polygon[:, i-1] )
                p = intersectionPointFromRays( yOne, ray )
                truncated = [truncated p]
            end
            # if the vertex is above y=1 and the next vertex is below y=1, find the intersection point and add it to the truncated polygon
            if polygon[2, i+1] < 1 # we know that i<size(polygon, 2), since the last vertex is on the x axis
                ray = Ray( polygon[:, i], polygon[:, i+1] - polygon[:, i] )
                q = intersectionPointFromRays( yOne, ray )
                truncated = [truncated q]
            end
        end
    end
    return truncated
end

"""
    subtract_piecewise(outerDome::Array{<:Rational,2}, innerDomes::Array{<:Array{<:Rational,2},1}) :: Array{<:Rational,2}

Compute the difference between a piecewise linear function `outerDome` and an array of piecewise linear functions `innerDomes`.
Each function is represented as a 2×n matrix of vertices in R^2. The first and last vertex of each function lie on the x-axis.
It is guaranteed that the `innerDomes` have disjoint support and are ordered along the x-axis.
The difference is only computed where the `innerDomes` are greater than zero.

# Arguments
- `outerDome::Array{<:Rational,2}`: The outer dome function, represented as a 2×n matrix of vertices.
- `innerDomes::Array{<:Array{<:Rational,2},1}`: The inner dome functions, represented as an array of 2×n matrices of vertices.

# Returns
- `diff_vertices::Array{Rational,2}`: The vertices of the difference function.
"""
function subtractPiecewise(outerDome::Array{<:Rational,2}, innerDomes::Array{<:Array{<:Rational,2},1}) :: Array{<:Rational,2}
    @assert !isempty( outerDome ) # outerDome should never be empty since we check for this in getCorridorIntersectionAreas
    # If there are no innerDomes return the outerDome
    if isempty( innerDomes ) return outerDome end
    # Initialize an empty matrix to store the vertices of the difference function
    diff_vertices = Matrix{Rational}(undef, 2, 0) 

    # Initialize index for outerDome
    outerDomeIndex = 1

    # Iterate over the array of innerDomes
    for innerDome in innerDomes
        # Initialize index for innerDome
        innerDomeIndex = 1

        # Find the index of the first vertex in outerDome that is not less than the first vertex in innerDome
        local nextOuterDomeIndex = findfirst( >=(innerDome[1, 1]), outerDome[1, :] )
        # If there is no such vertex, we can skip this innerDome
        if nextOuterDomeIndex === nothing  continue  end

        # Initialize diff_vertices_cache as an empty array
        diff_vertices_cache = SVector{2,Rational}[]

        # Handle vertices of outerDome that are before the current innerDome
        if nextOuterDomeIndex > outerDomeIndex
            diff_vertices = hcat( diff_vertices, outerDome[:, outerDomeIndex:nextOuterDomeIndex-1] )
            outerDomeIndex = nextOuterDomeIndex
        end

        while outerDomeIndex <= size( outerDome, 2 ) && innerDomeIndex <= size( innerDome, 2 )
            x_outer, y_outer = outerDome[:, outerDomeIndex]
            x_inner, y_inner = innerDome[:, innerDomeIndex]
        
            local vertex :: SVector{2,Rational}
            if x_outer < x_inner
                # If the x-value of outerDome is less than the x-value of innerDome, interpolate the y-value of innerDome
                y_at_x_outer = interp( innerDome, x_outer, innerDomeIndex )
                vertex = SA[ x_outer, safely( -, y_outer, max( 0, y_at_x_outer ) ) ]
                outerDomeIndex += 1
            elseif x_inner < x_outer
                # If the x-value of innerDome is less than the x-value of outerDome, interpolate the y-value of outerDome
                y_at_x_inner = interp( outerDome, x_inner, outerDomeIndex )
                vertex = SA[ x_inner, safely( -, y_at_x_inner, max( 0, y_inner ) ) ]
                innerDomeIndex += 1
            else
                # Otherwise we compute the difference directly
                vertex = SA[ x_outer, safely( -, y_outer, max( 0, y_inner ) ) ]
                outerDomeIndex += 1
                innerDomeIndex += 1
            end
            push!( diff_vertices_cache, vertex )

            if length( diff_vertices_cache )>=2
                local lastYDiff, secondLastYDiff = diff_vertices_cache[end][2], diff_vertices_cache[end-1][2]
                if (lastYDiff<0<secondLastYDiff  || lastYDiff>0>secondLastYDiff)
                    # an intersection has happened, so we need to replace the last entry with the intersection point
                    local direction = safely( -, diff_vertices_cache[end], diff_vertices_cache[end-1] )
                    local diffRay = Ray( diff_vertices_cache[end-1], direction )
                    local last = pop!( diff_vertices_cache )
                    local intersectionPoint = SA[getIntersectionTsWithXAxis( diffRay )[1], 0]
                    push!( diff_vertices_cache, intersectionPoint )
                    if diff_vertices_cache[end][2]>=0
                        push!( diff_vertices_cache, last )
                    end
                end
            end
        end
        
        diff_vertices = hcat( diff_vertices, reduce( hcat, diff_vertices_cache ) )
    end

    # If diff_vertices is empty, return the outerDome
    if isempty( diff_vertices ) return  outerDome  end

    # Handle remaining vertices of outerDome
    diff_vertices = hcat( diff_vertices, outerDome[:, outerDomeIndex:end] )

    return diff_vertices
end

"""
    interp(vertices::Array{<:Rational,2}, x::Rational, j::Int) :: Rational

Interpolate the y-value at a given x-value for a piecewise linear function represented by a 2×n matrix of vertices. The vertices must be ordered along the x-axis.

# Arguments
- `vertices::Array{<:Rational,2}`: The vertices of the piecewise linear function.
- `x::Rational`: The x-value at which to interpolate the y-value.
- `j::Int`: The current index of the vertices.

# Returns
- `y::Rational`: The interpolated y-value.
"""
function interp( vertices::Array{<:Rational,2}, x::Rational, j::Integer ) :: Rational
    # Check if x is within the range of the vertices
    if x < vertices[1, 1] || x > vertices[1, end]
        return 0
    end

    # Find the bracketing vertices
    x1, y1 = vertices[:, j-1]
    x2, y2 = vertices[:, j]

    # Interpolate the y-value
    y = try
        y1 + (x - x1) * (y2 - y1) // (x2 - x1)
    catch e
        if isa( e, OverflowError )
            big(y1) + (big(x) - big(x1)) * (big(y2) - big(y1)) // (big(x2) - big(x1))
        else
            rethrow( e )
        end
    end

    return y
end

"""
    compute_area(vertices::AbstractArray{<:Rational,2}) :: Rational{BigInt}

Compute the area under a piecewise linear function represented by a 2×n matrix of vertices. The vertices must be ordered along the x-axis.

# Arguments
- `vertices::AbstractArray{<:Rational,2}`: The vertices of the piecewise linear function.

# Returns
- `area::Rational{BigInt}`: The area under the function.
"""
function compute_area( vertices::AbstractArray{<:Rational,2} ) :: Rational{BigInt}
    # compute widths and heights
    ws = [ big(xNext) - big(x) for (x, xNext) in neighbors( vertices[1,:] ) ]
    hs = [ big(y) + big(yNext) for (y, yNext) in neighbors( vertices[2,:] ) ]

    # apply trapezoidal rule
    return ws ⋅ hs / 2
end

function find_intersections(dome1::Matrix{<:Rational}, dome2::Matrix{<:Rational})
    intersections = []

    rays1 = [Ray( p1, safely( -, p2, p1 ) ) for (p1, p2) in neighbors( eachcol( dome1 ) )]
    rays2 = [Ray( q1, safely( -, q2, q1 ) ) for (q1, q2) in neighbors( eachcol( dome2 ) )]

    for r1 in rays1, r2 in rays2
        # check if the x intervals overlap
        if r1.origin[1] > r2.origin[1] + r2.direction[1] || r2.origin[1] > r1.origin[1] + r1.direction[1]
            continue
        end
        # check if rays are parallel using determinant
        if r1.direction[1] * r2.direction[2] == r1.direction[2] * r2.direction[1]
            continue
        end

        t = getIntersectionTs( r1, r2 )

        if all( 0 .<= t .<= 1 )
            intersection_point = stepRayByT( r1, t[1] )
            push!( intersections, intersection_point )
        end
    end

    return intersections
end

"""
    alignDomes( dome1::Matrix{<:Rational}, dome2::Matrix{<:Rational} ) :: Matrix{Union{<:Rational, Tuple{<:Rational,<:Rational}}}

Align two piecewise linear functions represented by 2×n matrices of vertices along the x-axis.
The vertices must be ordered along the x-axis.
The domes must be non-empty.
"""
function alignDomes( dome1::Matrix{<:Rational}, dome2::Matrix{<:Rational} ) :: Matrix{Union{<:Rational, Tuple{<:Rational,<:Rational}}}
    intersections = find_intersections( dome1, dome2 )

    # initialize dictionary with intersections
    newDomeDict = Dict{Rational, Tuple{Rational,Rational}}()
    for (x, y) in intersections
        newDomeDict[x] = (y, y)
    end

    # add missing x-values to the dictionary
    for (x, y) in eachcol( dome1 )
        if !haskey( newDomeDict, x )
            # interpolate the y-value at x for dome2 and take the minimum of the two y-values
            i = findfirst( ≥(x), dome2[1, :] )
            if i === nothing || i == 1
                newDomeDict[x] = (y, zero( y ))
            else
                newDomeDict[x] = (y, interp( dome2, x, i ))
            end
        end
    end
    for (x, y) in eachcol( dome2 )
        if !haskey( newDomeDict, x )
            # interpolate the y-value at x for dome1 and take the minimum of the two y-values
            i = findfirst( ≥(x), dome1[1, :] )
            if i === nothing || i == 1
                newDomeDict[x] = (zero( y ), y)
            else
                newDomeDict[x] = (interp( dome1, x, i ), y)
            end
        end
    end

    newDome = reduce( hcat, [ [x, ys] for (x, ys) in newDomeDict ] )
    perm = sortperm( newDome[1, :] )
    return newDome[:, perm]
end

function removeRedundantVertices( dome::Matrix{<:Rational} ) :: Matrix{<:Rational}
    # remove redundant vertices with y=0, i.e. keep only indices with y=0 that
    # - are directly followed by a vertex with y>0
    # - are directly preceded by a vertex with y>0
    positive = dome[2, :] .> 0
    nextPositive = [positive[2:end]; false]
    prevPositive = [false; positive[1:end-1]]

    # Determine vertices to keep: those with y > 0 or adjacent to a vertex with y > 0
    keepVertex = positive .| nextPositive .| prevPositive

    # Apply the filter to the dome matrix
    dome = dome[:, keepVertex]

    return dome
end

function minDome( dome1::Matrix{<:Rational}, dome2::Matrix{<:Rational} ) :: Matrix{<:Rational}
    if isempty( dome1 ) || isempty( dome2 )
        return Matrix{Rational}(undef, 2, 0)
    end

    alignment = alignDomes( dome1, dome2 )
    dome = reduce( hcat, [ [x, minimum( ys )] for (x, ys) in eachcol( alignment ) ] )

    return removeRedundantVertices( dome )
end

function maxDome( dome1::Matrix{<:Rational}, dome2::Matrix{<:Rational} ) :: Matrix{<:Rational}
    if isempty( dome1 )
        return dome2
    elseif isempty( dome2 )
        return dome1
    end

    alignment = alignDomes( dome1, dome2 )
    dome = reduce( hcat, [ [x, maximum( ys )] for (x, ys) in eachcol( alignment ) ] )

    return removeRedundantVertices( dome )
end
