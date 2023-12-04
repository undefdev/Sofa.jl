using Sofa

using Test


const u = Point2d(1, 2).//1
const α = pi / 3
triple = (56, 33, 65)
sinα = triple[1] // triple[3]
cosα = triple[2] // triple[3]
tanα = sinα // cosα

# Calculate intersections
# const intersection_right = u[1]*cos(α) + u[1]*tan(α)*sin(α)
# const intersection_left  = u[2]*(-1)/tan(α)*cos(α) - u[2]*sin(α)
const intersection_right = u[1]*cosα + u[1]*tanα*sinα
const intersection_left  = u[2]*(-1)/tanα*cosα - u[2]*sinα

# Define two LRay that will intersect
xAxisRay = LRay( Point2d(0, 0).//1, Point2d(1, 0).//1 )
lray2 = rotate( LRay( u, Point2d(1, 0).//1 ), triple)

# Create array of LRay
lRays = [xAxisRay, lray2]
ts = getIntersectionTs(@view lRays[:])

lRays2 = [LRay(Point2d(-1, 2).//1, Point2d(1, 1//2).//1), LRay( Point2d(1,2).//1, Point2d(1, 1).//1 )]
ts2 = getIntersectionTs( @view lRays2[:] )

@testset "intersection ts" begin
    @test getIntersectionTs( Sofa.Ray( Point2d(0, 0).//1, Point2d(1, 0).//1 ), Sofa.Ray( Point2d(0, -1//2).//1, Point2d(0, 1).//1 ) ) == [0, 1//2].//1

    # Verify the intersection times for the rays
    # Check if the returned t-values are close to expected t-values, allowing a small error
    # The exact values you would compare here would depend on what you expect them to be
    @test ts[1, 2, 1, 1] == intersection_left
    @test ts[1, 2, 1, 2] == intersection_right

    local lrs = [LRay(Point2d(0,0).//1, Point2d(1,0).//1), LRay( Point2d(1//2,0).//1, Point2d(1, 2).//1 ), LRay( Point2d(1//4, 0).//1, Point2d(1,1).//1 )]
    @test getIntersectionTs( @view lrs[:] )[1,1:3,1,1] == [0, 1//2, 1//4]
end

@testset "testing truncation of polygon" begin
    @test Sofa.truncate( [ 0 2 3 4
                           0 2 3 0].//1 ) == [ 0 1 3+2//3 4
                                               0 1 1      0 ].//1 
    @test Sofa.truncate( [ -1  0    1    2   3 ;
                            0  2    0    2   0 ].//1 )  ==
                         [ -1 -1//2 1//2 1 3//2 5//2 3 ;
                            0  1    1    0 1    1    0 ].//1
end

@testset "outer dome computation" begin
    local lrs = @view lRays[:]
    local lrs2 = @view lRays2[:]
    @test outerDome( lrs, ts ) == (Matrix{Real}(undef, 0, 0), [], [])

    x, y, i, j = findOuterDomePeak( lrs2, ts2 )
    @test x == -1//3 && y == 2//3 && i==2 && j==1
    dome, domeIs, domeJs = outerDome( lrs2, ts2 )
    @test domeIs[]==i && domeJs[]==j
    @test dome[:,1]==[-1, 0] && dome[:,2]==[x, y] && dome[:,3]==[0, 0]

    dome = Rational{Int128}.(dome)
    local dome2 :: Matrix{Rational{Int128}} = [ -4 -1//2 -1//3 1 3 ;
                                                 0  2  2//3 1 0 ].//1

    display(dome)
    display(dome2)
    
    display(minDome( dome, dome2 ))
    display(maxDome( dome, dome2 ))
end

@testset "inner dome computation" begin
    local lrs2 = @view lRays2[:]
    components, indexComponents = innerDome( lrs2, ts2 )
    singleComponent = components[]
    @test singleComponent == [ -5 -1 -1//3 1 3 ;
                                0  2  2//3 2 0 ]
    @test innerDomeTruncated( lrs2, ts2 )[1][]== [ -5 -3 -1//2 -1//3 0 2 3 ;
                                                    0  1  1     2//3 1 1 0 ]
end

@testset "compute corridor intersection area" begin
    corridorSets = [
        [
            Corridor(Point2d( 0//1, 117//112 ), (33, 56, 65)),
            Corridor(Point2d( 2437//7280, 129//715 ), (56, 33, 65)),
            Corridor(Point2d( 49849//100555, 172919//101400 ), (119, 120, 169))
        ],
        [
            Corridor(Point2d( 0//1, 13//28 ), (33, 56, 65)),
            Corridor(Point2d( 1641//3640, 2793//1430 ), (56, 33, 65)),
            Corridor(Point2d( 33469//100555, 271199//101400 ), (119, 120, 169))
        ],
        [ # should have area 0
            Corridor(Point2d( 0//1, 897//1120 ), (33, 56, 65)),
            Corridor(Point2d( -156539//100555, -46573//101400 ), (119, 120, 169)),
            Corridor(Point2d( -267//7280, 4661//4290 ), (56, 33, 65)),
        ],
        [ # outer dome is inside inner dome, should have area 0 
            Corridor(Point2d( 0//1, 39//140 ), (33, 56, 65)),
            Corridor(Point2d( -235163//100555, 254819//101400 ), (119, 120, 169)),
            Corridor(Point2d( -22//65, 619//195 ), (56, 33, 65))
        ],
        [
            Corridor(Point2d( 0//1, 13//40 ), (33, 56, 65)),
            Corridor(Point2d( 7349//20111, 84467//101400 ), (119, 120, 169)),
            Corridor(Point2d( -7027//7280, 2753//2145 ), (56, 33, 65))
        ],
        [
            Corridor(Point2d( 0//1, 429//1120 ), (33, 56, 65)),
            Corridor(Point2d( -51707//100555, 248267//101400 ), (119, 120, 169)),
            Corridor(Point2d( -1461//1820, 7027//4290 ), (56, 33, 65))
        ],
        [
            Corridor(Point2d( 0//1, 767//1120 ), (33, 56, 65)),
            Corridor(Point2d( -6703//5915, 222059//101400 ), (119, 120, 169)),
            Corridor(Point2d( -495//364, 1779//1430 ), (56, 33, 65))
        ],
        [
            Corridor(Point2d( 0//1, 1209//1120 ), (33, 56, 65)),
            Corridor(Point2d( 46573//100555, 84467//101400 ), (119, 120, 169)),
            Corridor(Point2d( -257//260, 1957//4290 ), (56, 33, 65))
        ]
    ]
    
    for cs in corridorSets
        areas = getCorridorIntersectionAreas( cs )
        totalArea = sum( areas )
        outerLRays = lRaysFromCorridors( cs )[:,2]
        outerView = @view outerLRays[:]
        xAxisIntersectionTs = getIntersectionTsWithXAxis( outerView )
        leftX = maximum( xAxisIntersectionTs[:,1,1] )
        rightX = minimum( filter( x->x>=leftX, xAxisIntersectionTs[:,1,2] ) )
        sampledArea = sampleTotalCorridorIntersectionArea( cs, 1_000_000, leftX, rightX ) 
        @test isapprox( float(totalArea), sampledArea, atol=3e-2 )
    end

    triples = [
        (  7,   24,  25 ),
        (  33,  56,  65 ),
        ( 119, 120, 169 ),
        (  56,  33,  65 ),
        (  24,   7,  25 )
    ]
    for i in 1:1
    corridors = sampleCorridorsFromTriples( triples )
    display( "### sampled corridors" )
    display( corridors )
    display( "### computed corridors" )
    areas = getCorridorIntersectionAreas( corridors )
    display( "### computed areas" )
    display( areas )
    display( "### computed total area" )
    totalArea = sum( areas )
    display( totalArea )
    display( float(totalArea) )
    display( "### sampled total area")
    outerLRays = @view lRaysFromCorridors( corridors )[:,2]
    xAxisIntersectionTs = getIntersectionTsWithXAxis( outerLRays )
    leftX = maximum( xAxisIntersectionTs[:,1,1] )
    rightX = minimum( filter( x->x>=leftX, xAxisIntersectionTs[:,1,2] ) )
    sampledArea = sampleTotalCorridorIntersectionArea( corridors, 1_000_000, leftX, rightX ) 
    display(sampledArea)
    display( "deviation: $(float(totalArea) - sampledArea)")
    @test isapprox( float(totalArea), sampledArea, atol=1e-2 )
    end
end

@testset "initial subbox speedtest" begin
    triples = [
        ( 697, 696, 985 ),
        (  7,   24,  25 ),
        (  33,  56,  65 ),
        (  56,  33,  65 ),
        (  24,   7,  25 )
    ]    
    box = subBoxFromTriples( triples )
    sBox, upperBound = refineBounds(box, triples, 2.2195, 300, 2)
end

@testset "pythagorean triple generation" begin
    @test all( isPythagoreanTriple, genTriples() )
end

@testset "aligned refinement" begin
    triples = [
        ( 697, 696, 985 ),
        (  7,   24,  25 ),
        (  33,  56,  65 ),
        (  56,  33,  65 ),
        (  24,   7,  25 )
    ]    
    
    iTriples = intermediateTriples( triples )[[4]]
    insert!(iTriples, 1, triples[1] )
    box = subBoxFromTriples( triples )
 
    sBox = refineBounds( box, triples, 2.2195, 20, 2 )

    rBox = Sofa.refinementBoxesFromTriples( iTriples, 4, 1, 2.2195 )
    println("refining box with $(length(rBox)) boxes..." )
    iTriples = iTriples[2:end]
    Sofa.refineBoundsWithBoxes( sBox, triples, rBox, iTriples, 2.2195, 10, 2 )
end
