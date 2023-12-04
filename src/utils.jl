
"""
    neighbors(xs::Vector)

Return a vector of pairs of adjacent elements of `xs`. For example,
`neighbors([1,2,3])` returns `[(1,2), (2,3)]`.
"""
neighbors( xs::AbstractArray ) = zip(xs[1:end-1], xs[2:end])


function combineElements( a::Vector, b::Vector ) :: Vector{Vector}
    n, m = length( a ), length( b )
    ret = []
    for i in 0:m^n-1
        d = digits( i, base=m, pad=n ).+1
        push!( ret, collect( zip( a, b[d] ) ) )
    end
    return ret
end

safely( f::Function, x::Rational ) :: Rational = try
    f(x)
catch e
    if isa( e, OverflowError ) || isa( e, InexactError )
        f( big(x) )
    else
        rethrow( e )
    end
end

safely( op::Function, x::Rational, y::Rational ) :: Rational = try
    op(x,y)
catch e
    if isa( e, OverflowError ) || isa( e, InexactError )
        op( big(x), big(y) )
    else
        rethrow( e )
    end
end

safely( op::Function, x::SVector{2,<:Rational}, y::SVector{2,<:Rational} ) :: SVector{2,Rational} = try
    op(x,y)
catch e
    if isa( e, OverflowError ) || isa( e, InexactError )
        op( big.(x), big.(y) )
    else
        rethrow( e )
    end
end

safely( op::Function, x::SMatrix{2,2,<:Rational}, y::SVector{2,<:Rational} ) :: SVector{2,Rational} = try
    op(x,y)
catch e
    if isa( e, OverflowError ) || isa( e, InexactError )
        op( big.(x), big.(y) )
    else
        rethrow( e )
    end
end

safely( op::Function, x::SubArray{<:Rational, 1}, y::SubArray{<:Rational, 1} ) :: SVector{2,Rational} = try
    op(x,y)
catch e
    if isa( e, OverflowError ) || isa( e, InexactError )
        op( big.(x), big.(y) )
    else
        rethrow( e )
    end
end
