"""
Generate a regular lattice of points
$(TYPEDSIGNATURES)

# Arguments
- `dx` → element spacing on the x axis
- `dy` → element spacing on the y axis
- `ds` → displacement along x between rows of elements
- `f_cond::Function` → function of two arguments (`f(x,y) = ...`) returning `true` if element at position `x,y` must be kept and `false` otherwise

# Keyord Arguments
- `dx0 = 0` → x coordinate of the origin of the lattice
- `dy0 = 0` → y coordinate of the origin of the lattice
- `M::Int = 70` → Number of elements to generate per row of points before appliying the filtering function `f_cond`
- `N::Int = M` → Number of rows of points to generate before appliying the filtering function `f_cond`
"""
function generate_regular_lattice(dx::T, dy::T, ds::T, f_cond::Function = (x, y)->true;dx0 = T(0), dy0 = T(0), M::Int = 70,N::Int = M) where T<:Real
	# Function to generate x position as function of row,column number m,n
	x(m, n) = m * dx + n * ds + dx0
	# Function to generate y position as function of row,column number m,n
	y(n) = n * dy + dy0
	# Generate the elements. For each row, shift the columns to always have the search domain around x=0
	# out = [(x(m, n), y(n)) for n in -N:N for m in range(-M - round(n * ds / dx);length = 2 * M + 1) if f_cond(x(m, n), y(n))]
	out = [(x(m - round(Int,n * ds / dx), n), y(n)) for n in -N:N,m in -M:M if f_cond(x(m - round(Int,n * ds / dx), n), y(n))]
	return out
end
generate_regular_lattice(dx::Real,dy::Real,ds::Real,args...;kwargs...) = generate_regular_lattice(promote(dx,dy,ds)...,args...;kwargs...)

"""
    generate_rect_lattice(spacing_x::Real,spacing_y::Real[,f_cond];kwargs...)
# Summary
Generate a rectangular lattice of points (with different spacing among x and y directions)
# Arguments
- `spacing_x` → spacing between points on the x axis
- `spacing_y` → spacing between points on the y axis

See [`generate_regular_lattice`](@ref) for a description of `f_cond` and of  the keyword arguments
"""
generate_rect_lattice(spacing_x::Real,spacing_y::Real,args...;kwargs...) = generate_regular_lattice(spacing_x,spacing_y,0,args...;kwargs...)

"""
`generate_square_lattice(spacing::Real[,f_cond];kwargs...)`
# Summary
Generate a square lattice of points (with equal spacing among x and y directions)
# Arguments
- `spacing` → spacing between points on both x and y axis

See [`generate_regular_lattice`](@ref) for a description of `f_cond` and of  the keyword arguments
"""
generate_square_lattice(spacing::Real,args...;kwargs...) = generate_regular_lattice(spacing,spacing,0,args...;kwargs...)

"""
`generate_hex_lattice(spacing::Real[,f_cond];kwargs...)`
# Summary
Generate a hexagonal lattice of points (with equal distance between neighboring points).
The hexagonal lattice generated by this function has distance between points on the same
column √3 times greater than the distance between points on the same row.
# Arguments
- `spacing` → spacing between points

See [`generate_regular_lattice`](@ref) for a description of `f_cond` and of  the keyword arguments
"""
generate_hex_lattice(spacing::Real,args...;kwargs...) = generate_regular_lattice(spacing .* (1,√3/2,.5)...,args...;kwargs...)

# Get the conversion from linear to db and viceversa
"""
$(TYPEDSIGNATURES)
Convert a number from linear to dB
"""
lin2db(x::Real) = 10log10(x)
"""
$(TYPEDSIGNATURES)
Convert a number from dB to linear
"""
db2lin(x::Real) = 10^(x/10)

# Convert between frequency and wavelength
"""
$(TYPEDSIGNATURES)
Get the wavelength (in m) starting from the frequency (in Hz)
"""
f2λ(f::Real) = c₀/f
"""
$(TYPEDSIGNATURES)
Get the frequency (in Hz) starting from the wavelength (in m) 
"""
λ2f(λ::Real) = c₀/λ

## Here we have the functions for the coloring computation
function compute_F_cell(max_colours::Int;grid_max::Int=25)
    #=
    This function is used to compute all the possible 2x2 lattice generating matrices for possible coloring schemes up to 'max_colours' colors
    Computation is done with a brute-force approach, generating all possible 2x1 vectors with maximum elements up to grid_max
    =#
    
    # Find the current length of the pre-computed vector
    current_length = length(F_reuse_matrix.square)
    if current_length >= max_colours
        # We already computed the function for the required number of colors
        return
    end
    n_missing = max_colours - current_length
    append!(F_reuse_matrix.square,Vector{SMatrix{2,2,Float64,4}}(undef,n_missing))
    append!(F_reuse_matrix.triangular,Vector{SMatrix{2,2,Float64,4}}(undef,n_missing))
    # Compute the matrix for the square lattice
    _coloring_inner_bruteforce!(F_reuse_matrix.square,I,grid_max)
    # Compute the matrix for the triangular lattice
    _coloring_inner_bruteforce!(F_reuse_matrix.triangular,@SMatrix([1 0;cosd(60) sind(60)]),grid_max)
end

function _coloring_inner_bruteforce!(T_mat, rot_mat, grid_max)
    max_colours = length(T_mat)
    check_vec = fill((typemax(Int),typemax(Int)),max_colours)
    @inline norm2(x) = sum(abs2.(x))
    @inbounds for x1 = 0:grid_max, y1 = 0:grid_max, x2 = -grid_max:grid_max, y2 = 0:grid_max
        mat = @SMatrix [x1 y1;x2 y2]
        # Compute the determinant
        t_det = Int(det(mat))
        # Skip points which are not likely to give useful results
        if (t_det < 1)  || (t_det > max_colours) || (t_det > grid_max^2) || (maximum(abs.(mat)) > ceil(sqrt(t_det) + 3))
            continue
        end
        # Compute the angle between the basis vectors identified by [x1,y1] and [x2,y2]
        angle = abs(atan(y1,x1) - atan(y2,x2))
        # Skip cases where the angle between vectors is either too acute or too obtuse
        if abs(π/2-angle) > π/4
            continue
        end
        # Create temp variables for computation of the minimum distance
        dmat = mat*rot_mat
        # Compute frobenius norm and minimum squared distance for the candidate lattice generating matrix
        frobe = round(Int,norm2(dmat))
        # display(frobe)
        # Minimum squared distance is either the modulo of one of the vectors or the modulo of sum or difference of them
        dmin = round(Int,minimum((norm2(dmat[1,:]), norm2(dmat[2,:]), norm2(sum(dmat;dims=1)), norm2(diff(dmat;dims=1)))))
        # Check if the current realization is better than the saved one
        if isless((-dmin,frobe),check_vec[t_det])
            # Update the check_vec
            check_vec[t_det] = (-dmin,frobe)
            # Update the vector containing the generating matrices
            T_mat[t_det] = mat'
        end
    end
end

function generate_F_reuse_matrix(lattice_type::Symbol=:triangular,N_Colours::Int=4;max_colours::Int=max(10,N_Colours))
    compute_F_cell(max_colours)
    return getproperty(F_reuse_matrix,lattice_type)[N_Colours]
end

"""
    generate_colors(BeamCenters::AbstractVector,N_Colours::Int=4;lattice_type::Symbol=:triangular)   

Provide a the coloring breakdown for a given set of lattice points.
# Arguments
- `BeamCenters` → Vector of Tuple or StaticVectors expressing the U-V coordinates of each point in the lattice for which the coloring is being computed
- `N_Colours` → Number of colors to divide the lattice in. Defaults to `4`

# keyword Arguments
- `lattice_type` → Symbol that can either be `:triangular` or `:square`, idenifying the type of lattice. Defaults to `:triangular`
"""
function generate_colors(BeamCenters::AbstractVector,N_Colours::Int=4;first_color_coord=nothing,first_color_idx=nothing,precision_digits::Int=7,lattice_type::Symbol=:triangular)
    #=
    **************************************************************************
       Generate frequency colouring file
    **************************************************************************

     References:
     [1]   "On the frequency allocation for mobile radio telephone systems", C. de
           Almeida; R. Palazzo, Proceedings of 6th International Symposium on
           Personal, Indoor and Mobile Radio Communications, Year: 1995, Volume: 1, Pages: 96 - 99 vol.1
     [2]   P. Angeletti, "Simple implementation of vectorial modulo operation based
           on fundamental parallelepiped," in Electronics Letters, vol. 48, no. 3, pp. 159-160,
           February 2 2012. doi: 10.1049/el.2011.3667

     Authors: Alberto Mengali, 2020, European Space Agency

    Input Arguments:
    BeamCenters:         A vector containing the U,V coordinates of the beam centers as tuples or staticvectors
    N_Colours:           An integer number specifying the number of colors to generate in the association
    first_color_coord:   A tuple or static vector containing the U,V coordinates of the beam that will have the first color
    first_color_idx:     The beam index of the beam containing the first color, either this variable of first_order_coord can be specified, not together
    precision_digits:    The number of digits to be used in the flooring operation
    =#

    # Check if either the first color coordinates or first color idx are given
    if first_color_coord === nothing
        if first_color_idx === nothing
            # Initialize to coorinates of the first beam
            first_color_idx = 1;
        end
        first_color_coord = BeamCenters[first_color_idx]
    else
        if first_color_idx !== nothing
            @warn "Both first_color_idx and first_color_coord were given, disregarding the idx variable"
        end
    end
    if lattice_type ∉ (:triangular, :square)
        @error "The specified lattice type ($lattice_type) is not recognized, it should be either :triangular or :square"
    end

    n_beams = length(BeamCenters)

    # If only 1 color is requested, simply return the trivial result
    if N_Colours == 1
        Colours = ones(n_beams)
        nbeams_per_color = n_beams
        idxs = fill(fill(true,nbeams))
    end

    # Find the minimum distance in U and V
    minU_dist = minimum(diff(sort(first.(BeamCenters)) |> x -> unique!(y -> round(y;digits=precision_digits),x)))
    minV_dist = minimum(diff(sort(last.(BeamCenters)) |> x -> unique!(y -> round(y;digits=precision_digits),x)))
    if lattice_type === :triangular
        beamU_spacing = 2minU_dist
        beamV_spacing = 2minV_dist
        # Matrix to normalize the grid points in u-v into a integer grid with y(v) axis not pointing north but north-east with 60° inclination
        D = @SMatrix [1 -1/2;0 1]
    elseif lattice_type === :square
        beamU_spacing = minU_dist
        beamV_spacing = minV_dist
        D = @SMatrix [1 0;0 1]
    end
    # Get the coloring generating matrix
    F_reuse_matrix = generate_F_reuse_matrix(lattice_type,N_Colours)
    # Create the set that will contain the unique results
    unique_colors_vecs = SVector{2,Int}[]
    # Initialize the colors vector
    Colors = similar(BeamCenters, Int)
    @inbounds for (n,p₀) ∈ enumerate(BeamCenters)
        # Compute the integer beam centers
        p₀_normalized = p₀ ./ SA_F64[beamU_spacing, minV_dist]
        # Transform the beam centers from u-v coordinates in radians into integer indexes
        beam_index_vector = round.(Int,D*p₀_normalized);
        # Find the values of the beam indexes modulo F_reuse_matrix
        unique_beam_index = round.(Int,beam_index_vector .- (F_reuse_matrix*floor.(round.(inv(F_reuse_matrix)*beam_index_vector,digits=precision_digits))))
        # Check if this color has already been assigned/registered
        idx = findfirst(x -> x == unique_beam_index,unique_colors_vecs)
        if idx isa Nothing
            push!(unique_colors_vecs, unique_beam_index)
            cidx = length(unique_colors_vecs)
        else
            cidx = idx
        end
        Colors[n] = cidx
    end
    
    return Colors
end