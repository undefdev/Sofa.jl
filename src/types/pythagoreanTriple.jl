
PythagoreanTriple = NTuple{3,Int}
TripleArray = Array{PythagoreanTriple,1}

angleFromTriple( t::PythagoreanTriple ) = asin( t[1]//t[3] )

mirror( t::PythagoreanTriple ) :: PythagoreanTriple = (t[2], t[1], t[3])

isPythagoreanTriple( t::PythagoreanTriple ) = t[1]^2 + t[2]^2 == t[3]^2

isPrimitive( t::PythagoreanTriple ) = gcd( t... ) == 1

function genTriples( paramMax::Int=31, maxTriple::PythagoreanTriple=(84, 13, 85), cutoff=1000 ) :: TripleArray
    # generate all primitive triples (a, b, c) with m,n,k <= paramMax and a<b<c<cutoff
    rawSorted = [ tuple( sort!( [k*(m^2-n^2), k*2*m*n, k*(m^2+n^2)] )... ) for k in 1:paramMax for n in 1:paramMax for m in n+1:paramMax ]
    capped = filter( t -> t[3] < cutoff, rawSorted )
    deduplicated = unique( capped )
    primitive = filter( isPrimitive, deduplicated )
    # sort triples by angle from smallest to largest
    sort!( primitive, by=angleFromTriple )
    # mirror all triples so far up to maxTriple and append them to the list
    i = findfirst( t -> t==mirror(maxTriple), primitive )
    return vcat( primitive , reverse( mirror.( primitive[i:end] ) ) )
end

function intermediateTriple( t1::PythagoreanTriple, t2::PythagoreanTriple, refTriples::TripleArray ) :: PythagoreanTriple
    i1 = findfirst( t -> t==t1, refTriples )
    i2 = findfirst( t -> t==t2, refTriples )
    midAngle = ( angleFromTriple(t1) + angleFromTriple(t2) ) / 2
    # get index of triple with angle closest to midAngle
    iMid = findmin( i -> abs( angleFromTriple(refTriples[i]) - midAngle ), i1+1:i2-1 )[2] + i1
    return refTriples[iMid]
end

function intermediateTriples( triples::TripleArray ) :: TripleArray
    sortedTriples = sort( triples, by=angleFromTriple ) # first triple is middle triple by convention
    minTriple, maxTriple = sortedTriples[1], sortedTriples[end]
    referenceTriples = genTriples()
    extendedTriples = sortedTriples
    if minTriple != referenceTriples[1]
        insert!( extendedTriples, 1, referenceTriples[1] )
    end
    if maxTriple != referenceTriples[end]
        push!( extendedTriples, referenceTriples[end] )
    end
    return ( x -> intermediateTriple( x[1], x[2], referenceTriples )).( neighbors(extendedTriples) )
end
