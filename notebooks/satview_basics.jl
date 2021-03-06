### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 35da0c45-6d91-4ad9-811f-5d633684208e
begin
	using StaticArrays
	using LinearAlgebra
	using Unitful
	using Unitful.DefaultSymbols
	using SatelliteToolbox: Ellipsoid
	using Parameters
	import Proj4
	using DocStringExtensions
end

# ╔═╡ c3a5e97e-00ed-44af-9c67-fa8198900fbd
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	using PlutoDevMacros
	using PlutoUtils
	using PlutoTest
end
  ╠═╡ =#

# ╔═╡ e2f50ace-ad14-4e7c-af74-abf3ca9df9fb
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Packages
"""
  ╠═╡ =#

# ╔═╡ 82c0d6f9-373d-4866-81eb-9cd2d0981310
# ╠═╡ skip_as_script = true
#=╠═╡
ToC()
  ╠═╡ =#

# ╔═╡ 2e6bcf9b-002f-4fbc-80c1-0cfd2ab1d253
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Exports
"""
  ╠═╡ =#

# ╔═╡ 1e3da0c9-2f96-4637-9093-ac7f10c1ad27
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Helper Functions
"""
  ╠═╡ =#

# ╔═╡ 48c73104-fe4c-4543-942a-6f23b0fd2547
const Point{N} = Union{Tuple{Vararg{<:Number, N}},StaticVector{N,<:Number}}

# ╔═╡ 48c73104-fe4c-4543-0001-6f23b0fd2547
const Point2D = Point{2}

# ╔═╡ 48c73104-fe4c-4543-0002-6f23b0fd2547
const Point3D = Point{3}

# ╔═╡ c4402f72-67ac-4630-a651-da81c1df71bf
"""
	ExtraOutput
Struct used inside `SatView` to use multiple dispatch to give more than the standard output. 
See for example [`get_range`](@ref)
"""
struct ExtraOutput end

# ╔═╡ e796d5be-e011-45d3-ad90-58769feb5e85
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Show/Print
"""
  ╠═╡ =#

# ╔═╡ d890aff1-dbd0-451c-bf14-bde9758c3be0
function _print_angle(io,val,displayname,last=false)
	print(io,"$displayname=")
	print(io,round(val;digits=2) * rad)
	print(io," (")
	print(io,round(rad2deg(val);digits=2) * °)
	print(io,")")
	last || print(io,", ")
end

# ╔═╡ fdbbc8d9-83d6-413e-aff8-aca18f24dfea
function _print_length(io,val,displayname,last=false)
	print(io,"$displayname=")
	mval = val < 1000 ? round(val;digits=2) * m : round(val/1000;digits=2) * km
	print(io,mval)
	last || print(io,", ")
end

# ╔═╡ cc0cae75-ba10-4a62-b0ef-22259e40a083
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Spherical Ellipsoid
"""
  ╠═╡ =#

# ╔═╡ 855553e3-491b-4e8a-a482-95855697e063
"""
	SphericalEllipsoid(r = 6371e3)
Define a spherical ellipsoid of radius `r` to be used for the various transformation in [`SatView`](@ref). Defaults to 6371km radius
"""
SphericalEllipsoid(r = 6371e3) = Ellipsoid(r,0.0)

# ╔═╡ 1f7a7093-ce8e-461f-8b91-69266de86748
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## geod_geodesic
"""
  ╠═╡ =#

# ╔═╡ f83fc65f-5f7b-498d-9ed3-0762565ad710
Proj4.geod_geodesic(e::Ellipsoid) = Proj4.geod_geodesic(e.a, e.f)

# ╔═╡ 4714c6ae-27d9-47db-b12e-126283b10606
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# EarthModel
"""
  ╠═╡ =#

# ╔═╡ f9cc2880-bab1-4be3-b67b-d43508df8d3b
begin
"""
$(TYPEDEF)
Geometrical model of the Earth ellipsoid to be used for all the view angles between satellite and points on earth.
# Fields
$(TYPEDFIELDS)

A single instance of this structure should be used for all satellites of a given simulation/constellation.\\
Changes to any of the two fields (via `setproperty!`) will trigger an automatic recomputation of the other field.\\
When called without arguments, it defaults to a spherical earth with a radius of 6371 km.

See also: [`Ellipsoid`](@ref), [`SphericalEllipsoid`](@ref), [`SatView`](@ref)
"""
@with_kw mutable struct EarthModel
	"Ellipsoid structure, used for the various point of view conversion."
	ellipsoid::Ellipsoid{Float64} = SphericalEllipsoid()
	"Extended geod structure, used for the inverse geodesic problem to compute distance and azimuth between points on earth. (Relies on GeographicLib)." 
	geod::Proj4.geod_geodesic = Proj4.geod_geodesic(ellipsoid)
end
EarthModel(e::Ellipsoid) = EarthModel(;ellipsoid = e)
end

# ╔═╡ 1f2b9b74-de46-401d-8a46-0434b9f9aca1
function Base.setproperty!(value::EarthModel, name::Symbol, x)
	if name === :ellipsoid
		setfield!(value, name, x)
		setfield!(value, :geod, Proj4.geod_geodesic(x))
	else
		setfield!(value, name, x)
		setfield!(value, :ellipsoid, Ellipsoid(x.a, x.f))
	end
end

# ╔═╡ 6a5cb372-60cb-4ffc-b4f0-22e4016104e7
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Angle Types
"""
  ╠═╡ =#

# ╔═╡ e832c1b7-8c04-4146-90f6-1628e91fea2a
const UnitfulAngleType = Union{typeof(°),typeof(rad)}

# ╔═╡ 64cf1b8b-6686-4eeb-a3cc-786300ea7c7d
const UnitfulAngleQuantity = Quantity{<:Real,<:Any,<:UnitfulAngleType}

# ╔═╡ 1d023a0c-a93a-451c-a894-1d1f6a4b78a9
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# SatViewCoordinate types
"""
  ╠═╡ =#

# ╔═╡ b8ce87d4-4768-4e8a-a16e-9e68b00b6617
abstract type SatViewCoordinate end

# ╔═╡ e3c221c6-6c4a-4b5f-93a6-30d3508ac9d2
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## LLA
"""
  ╠═╡ =#

# ╔═╡ 8368ae01-ce53-449e-87dd-8dfa3f29f8f4
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
Here we want to define a structure that contains useful informations and functions to perform conversions between the view from the satellite based on it's orbital position and points on ground
"""
  ╠═╡ =#

# ╔═╡ f951805e-515a-475f-893f-bb8b968e425c
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## ERA
"""
  ╠═╡ =#

# ╔═╡ 86ae20a9-e69c-4d63-9119-395449e9ac09
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
ERA stands for Elevation, Range and Azimuth and is used to express the position of a satellite relative to an observer in spherical coordinates.
The elevation is the angle of the pointing with respect to the local horizon of the observer, meaning the plane where the observer is located that is perpendicular to the gravity vector acting on the observe (or in an alternative definition, the plane where the observer is located that is parallel to the tangent plane to the earth ellipsoid at the given lat and lon positions of the observer.
"""
  ╠═╡ =#

# ╔═╡ 9be2fd5c-4b6c-4e13-b2aa-fb7120a504b7
begin
	"""
	Elevation, Range and Azimuth for a target point on space as seen from a source point on or above the earth surface
	
	# Fields
	- `el::Float64`: Elevation view angle (`0 <= el <= π/2`) between source and target point [rad].
	- `r::Float64`: Range (`r >= 0`) between the source and target points [m].
	- `az::Float64`: Azimuth view angle between source and target point [rad], computed from West to North from the source point perspective. Values provided are automatically converted between -π and π
	
	# Constructors
	
		ERA(el::Real,r::Real,az::Real)
		ERA(el::UnitfulAngleQuantity,r::Real,az::Real)
		ERA(el::UnitfulAngleQuantity,r::Real,az::UnitfulAngleQuantity)
		ERA(el,r::Unitful.Length,az)
	
	where `UnitfulAngleQuantity` is a `Unitful.Quantity` of unit either `u"rad"` or `u"°"`.
	"""
	@with_kw_noshow struct ERA <: SatViewCoordinate
		el::Float64 # Elevation in radians
		r::Float64 # Range in meters
		az::Float64 # Azimuth in radians
		
		function ERA(el::Real,r::Real,az::Real)
			(isnan(el) || isnan(r) || isnan(az)) && return new(NaN,NaN,NaN)  
			@assert el >= 0 && el <= π/2 "Elevation (and azimuth) should be provided in radians and between 0 and π/2.  If you wish to use degrees as inputs please multiply the first and last arguments by `°` (from Unitful.jl), which is re-exported by TelecomUtils and can be written with `\degree<TAB>`"
			@assert r >= 0 "Range must be positive"
			new(el,r,rem2pi(az,RoundNearest))
		end
	end

	# Define a constructor that takes combinations of real numbers and angles/lengths
	ERA(el::UnitfulAngleQuantity,r::Real,az::Real) = ERA(
		uconvert(u"rad",el) |> ustrip,
		r,
		az,
	)
	ERA(el::UnitfulAngleQuantity,r::Real,az::UnitfulAngleQuantity) = ERA(
		el,
		r,
		uconvert(u"rad",az) |> ustrip,
	)
	ERA(el,r::Unitful.Length,az) = ERA(
		el,
		uconvert(u"m",r) |> ustrip,
		az,
	)
	
	# Show
	function Base.show(io::IO,era::ERA)
		print(io,"ERA(")
		_print_angle(io,era.el,"el",false)
		_print_length(io,era.r,"r",false)
		_print_angle(io,era.az,"az",true)
		print(io,")")
	end
end

# ╔═╡ 16782c72-ecb1-48ec-8510-78e2e0689a10
function Base.isapprox(x1::ERA, x2::ERA)
	x1.el ≉ x2.el && return false
	# Don't care about different azimuth if elevation is 90°
	x2.el ≉ π/2 && x1.az ≉ x2.az && return false
	x1.r ≉ x2.r && return false
	return true
end

# ╔═╡ 7344190c-7989-4b55-b7be-357f7d6b7370
Base.isnan(era::ERA) = isnan(era.el)

# ╔═╡ f207d849-ebff-4e6c-95bb-50693cb7c9b6
begin
	"""
	Identify a point on or above earth using geodetic coordinates
	
	# Fields
	- `lat::Float64`: Latitude (`-π/2 <= lat <= π/2`) of the point [rad].
	- `lon::Float64`: Longitude of the point [rad].
	- `alt::Float64`: Altitude of the point above the reference earth ellipsoid [m].
	
	# Constructors
		LLA(lat::Real,lon::Real,alt::Real)
		LLA(lat::UnitfulAngleQuantity,lon::Real,alt::Real)
		LLA(lat::UnitfulAngleQuantity,lon::UnitfulAngleQuantity,alt::Real)
		LLA(lat,lon,alt::Unitful.Length)
		LLA(lat,lon) # Defaults to 0.0 altitude
	
	where `UnitfulAngleQuantity` is a `Unitful.Quantity` of unit either `u"rad"` or `u"°"`.
	"""
	@with_kw_noshow struct LLA <: SatViewCoordinate
		lat::Float64 # Latitude in radians
		lon::Float64 # Longitude in radians
		alt::Float64 # Altitude in meters
		
		function LLA(lat::Real,lon::Real,alt::Real)
			(isnan(lat) || isnan(lon) || isnan(alt)) && return new(NaN,NaN,NaN)  
			l2 = rem2pi(lon,RoundNearest)
			@assert abs(lat) <= π/2 "Latitude should be given in radians and between -π/2 and π/2. If you wish to use degrees as inputs please multiply the first two arguments by `°` (from Unitful.jl), which is re-exported by TelecomUtils and can be written with `\degree<TAB>`"
			new(lat,l2,alt)
		end
	end
	
	# Constructor without altitude, assume it is 0
	LLA(lat,lon) = LLA(lat,lon,0.0)
	
	# Define a constructor that takes combinations of real numbers and angles/lengths
	LLA(lat::UnitfulAngleQuantity,lon::Real,alt::Real) = LLA(
		uconvert(u"rad",lat) |> ustrip,
		lon,
		alt)
	LLA(lat::UnitfulAngleQuantity,lon::UnitfulAngleQuantity,alt::Real) = LLA(
		lat,
		uconvert(u"rad",lon) |> ustrip,
		alt)
	LLA(lat,lon,alt::Unitful.Length) = LLA(
		lat,
		lon,
		uconvert(u"m",alt) |> ustrip)	
end

# ╔═╡ 3033d5c1-d8e0-4d46-0001-7dec4ff7afbd
# Show method for LLA
function Base.show(io::IO,lla::LLA)
	print(io,"LLA(")
	_print_angle(io,lla.lat,"lat",false)
	_print_angle(io,lla.lon,"lon",false)
	_print_length(io,lla.alt,"alt",true)
	print(io,")")
end

# ╔═╡ 3033d5c1-d8e0-4d46-a4b9-7dec4ff7afbd
function Base.isapprox(x1::LLA, x2::LLA)
	x1.alt ≉ x2.alt && return false
	# Don't care about different longitude if latitude is ±90°
	abs(x1.lat) ≈ π/2 && abs(x2.lat) ≈ π/2 && return true
	# Return true if all the lat and lon are matching
	x1.lat ≈ x2.lat && (x1.lon ≈ x2.lon || abs(x1.lon) ≈ abs(x2.lon) ≈ π) && return true
	return false
end

# ╔═╡ 528beffe-0707-4661-8970-def1b1e00ea5
Base.isnan(lla::LLA) = isnan(lla.lat)

# ╔═╡ 11938cb6-46b3-0001-96c0-ef6424d1d0db
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Inverse Geodesic Problem
"""
  ╠═╡ =#

# ╔═╡ 11938cb6-46b3-0002-96c0-ef6424d1d0db
"""
    geod_inverse(geod::Proj4.geod_geodesic, lonlat1::AbstractVector{Cdouble}, lonlat2::AbstractVector{Cdouble})
	geod_inverse(geod::Proj4.geod_geodesic, lla1::LLA, lla2::LLA)
Solve the inverse geodesic problem.

# Args

- `g`       → the geod_geodesic object specifying the ellipsoid.
- `lonlat1` → point 1 (degrees), where lat ∈ [-90, 90], lon ∈ [-540, 540) 
- `lonlat2` → point 2 (degrees), where lat ∈ [-90, 90], lon ∈ [-540, 540) 

# Outputs

- `dist` → distance between point 1 and point 2 (meters).
- `azi1` → azimuth at point 1 (degrees) ∈ [-180, 180)
- `azi2` → (forward) azimuth at point 2 (degrees) ∈ [-180, 180)

# Remarks

If either point is at a pole, the azimuth is defined by keeping the longitude fixed,
writing lat = 90 +/- eps, and taking the limit as eps -> 0+.
"""
function geod_inverse(geod::Proj4.geod_geodesic, lonlat1::AbstractVector{Cdouble}, lonlat2::AbstractVector{Cdouble})
	dist = Ref{Cdouble}()
	azi1 = Ref{Cdouble}()
	azi2 = Ref{Cdouble}()
	ccall((:geod_inverse, Proj4.libproj), Cvoid, (Ptr{Cvoid},Cdouble,Cdouble,Cdouble,
			Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
			pointer_from_objref(geod), lonlat1[2], lonlat1[1], lonlat2[2], lonlat2[1], dist, azi1, azi2)
	dist[], azi1[], azi2[]
end

# ╔═╡ 11938cb6-46b3-0003-96c0-ef6424d1d0db
function geod_inverse(geod::Proj4.geod_geodesic, lla1::LLA, lla2::LLA)
	lonlat1 = rad2deg.(SA_F64[lla1.lon,lla1.lat])
	lonlat2 = rad2deg.(SA_F64[lla2.lon,lla2.lat])
	geod_inverse(geod,lonlat1,lonlat2)
end

# ╔═╡ 4c06d21c-ac14-4522-bf25-2e0a1ed2d6b9
begin
	export Ellipsoid, SphericalEllipsoid
	export EarthModel
	export UnitfulAngleQuantity, UnitfulAngleType, °
	export km
	export LLA, ERA
	export geod_inverse
end

# ╔═╡ 11938cb6-46b3-499b-96c0-ef6424d1d0db
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Tests
"""
  ╠═╡ =#

# ╔═╡ d2c248b1-c48e-437b-a910-edcc59b4424f
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## LLA
"""
  ╠═╡ =#

# ╔═╡ 7b306ed5-4bda-465d-abf2-4d07cb4642c1
# ╠═╡ skip_as_script = true
#=╠═╡
@test LLA(10°,10°,1000) ≈ LLA((10+100*eps())*°,10°,1000)
  ╠═╡ =#

# ╔═╡ c45f6aac-bff3-4c98-8bd5-c91a98c9eef7
# ╠═╡ skip_as_script = true
#=╠═╡
@test LLA(90°,10°,1000) ≈ LLA(90°,130°,1000)
  ╠═╡ =#

# ╔═╡ 1f1505a6-c6af-4fdd-9583-6e783e76de4f
# ╠═╡ skip_as_script = true
#=╠═╡
@test LLA(40°,-180°,1000) ≈ LLA(40°,180°,1000)
  ╠═╡ =#

# ╔═╡ 980e48bd-9dd1-4195-8086-40785a7f43e1
# ╠═╡ skip_as_script = true
#=╠═╡
@test LLA(10°,10°,1000) !== LLA((10+100*eps())*°,10°,1000)
  ╠═╡ =#

# ╔═╡ fc146750-46b2-4084-a6d2-0d91c0e104e6
# ╠═╡ skip_as_script = true
#=╠═╡
LLA(1,1,NaN) |> isnan
  ╠═╡ =#

# ╔═╡ fc3e96a6-7af8-4ef0-a4ca-b33509aa512f
md"""
## ERA
"""

# ╔═╡ ed2f2fa3-be19-4e2d-ad28-36cb053ed6bd
# ╠═╡ skip_as_script = true
#=╠═╡
@test ERA(10°,1000,20°) == ERA(10°,1km,deg2rad(20)*rad)
  ╠═╡ =#

# ╔═╡ 6151e4d1-a7a6-4fa7-bc77-0e7f3b2a3cc0
# ╠═╡ skip_as_script = true
#=╠═╡
@test ERA(90°,1000,20°) ≈ ERA(90°,1km,deg2rad(90)*rad)
  ╠═╡ =#

# ╔═╡ b5449c6a-9c43-4f9c-81f7-f01277acb109
# ╠═╡ skip_as_script = true
#=╠═╡
ERA(1,1,NaN) |> isnan
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DocStringExtensions = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
PlutoDevMacros = "a0499f29-c39b-4c5c-807c-88074221b949"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUtils = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
Proj4 = "9a7e659c-8ee8-5706-894e-f68f43bc57ea"
SatelliteToolbox = "6ac157d9-b43d-51bb-8fab-48bf53814f4a"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
DocStringExtensions = "~0.8.6"
Parameters = "~0.12.3"
PlutoDevMacros = "~0.4.5"
PlutoTest = "~0.2.2"
PlutoUtils = "~0.5.9"
Proj4 = "~0.7.6"
SatelliteToolbox = "~0.9.4"
StaticArrays = "~1.3.3"
Unitful = "~1.10.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0-rc1"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.Chain]]
git-tree-sha1 = "8c4920235f6c561e401dfe569beb8b924adad003"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.5.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9489214b993cd42d17f44c36e359bf6a7c919abf"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "924cdca592bc16f14d2f7006754a621735280b74"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.1.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "681ea870b918e7cff7111da58791d7f718067a19"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "9267e5f50b0e12fdfd5a2455534345c4cf2c7f7a"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.14.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Glob]]
git-tree-sha1 = "4df9f7e06108728ebf00a0a11edee4b29a482bb2"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.0"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b7bc05649af456efc75d178846f47006c2c4c3c7"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.6"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "ec2e30596282d722f018ae784b7f44f3b88065e4"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.6"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OptionalData]]
git-tree-sha1 = "d047cc114023e12292533bb822b45c23cb51d310"
uuid = "fbd9d27c-2d1c-5c1c-99f2-7497d746985d"
version = "1.0.0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PROJ_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "LibSSH2_jll", "Libdl", "Libtiff_jll", "MbedTLS_jll", "Pkg", "SQLite_jll", "Zlib_jll", "nghttp2_jll"]
git-tree-sha1 = "2435e91710d7f97f53ef7a4872bf1f948dc8e5f8"
uuid = "58948b4f-47e0-5654-a9ad-f609743f8632"
version = "700.202.100+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0044b23da09b5608b4ecacb4e5e6c6332f833a7e"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoDevMacros]]
deps = ["MacroTools", "Requires"]
git-tree-sha1 = "994167def8f46d3be21783a76705228430e29632"
uuid = "a0499f29-c39b-4c5c-807c-88074221b949"
version = "0.4.5"

[[deps.PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "17aa9b81106e661cffa1c4c36c17ee1c50a86eda"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8d1f54886b9037091edf146b517989fc4a09efec"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.39"

[[deps.PlutoUtils]]
deps = ["AbstractPlutoDingetjes", "Chain", "Colors", "DocStringExtensions", "Glob", "HypertextLiteral", "OrderedCollections", "PlutoDevMacros", "PlutoUI", "PrettyTables", "Reexport", "Requires", "StaticArrays", "UUIDs"]
git-tree-sha1 = "3f8dfe27dbb980ad5e83ecd641ded8eed91f3265"
uuid = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
version = "0.5.9"

[[deps.PolynomialRoots]]
git-tree-sha1 = "5f807b5345093487f733e520a1b7395ee9324825"
uuid = "3a141323-8675-5d76-9d11-e1df1406c778"
version = "1.0.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Proj4]]
deps = ["CEnum", "CoordinateTransformations", "PROJ_jll", "StaticArrays"]
git-tree-sha1 = "5f15f1c647b563e49f655fbbfd4e2ade24bd3c64"
uuid = "9a7e659c-8ee8-5706-894e-f68f43bc57ea"
version = "0.7.6"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.ReferenceFrameRotations]]
deps = ["Crayons", "LinearAlgebra", "Printf", "StaticArrays"]
git-tree-sha1 = "d526371cec370888f485756a4bf8284ab531860b"
uuid = "74f56ac7-18b3-5285-802d-d4bd4f104033"
version = "1.0.1"

[[deps.RemoteFiles]]
deps = ["Dates", "FileIO", "HTTP"]
git-tree-sha1 = "54527375d877a64c55190fb762d584f927d6d7c3"
uuid = "cbe49d4c-5af1-5b60-bb70-0a60aa018e1b"
version = "0.4.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SQLite_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "b6d006c4c57278d532de38912e16adf626c949c7"
uuid = "76ed43ae-9a5d-5a62-8c75-30186b810ce8"
version = "3.38.4+0"

[[deps.SatelliteToolbox]]
deps = ["Crayons", "Dates", "DelimitedFiles", "Interpolations", "LinearAlgebra", "OptionalData", "Parameters", "PolynomialRoots", "PrettyTables", "Printf", "Reexport", "ReferenceFrameRotations", "RemoteFiles", "SparseArrays", "StaticArrays", "Statistics"]
git-tree-sha1 = "1831cced8785398bf38577e8cf46380d349cf4c9"
uuid = "6ac157d9-b43d-51bb-8fab-48bf53814f4a"
version = "0.9.4"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "74eaf352c0cef1e32ce7394bcc359d9199a28cf7"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.3.6"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "b95e0b8a8d1b6a6c3e0b3ca393a7a285af47c264"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.10.1"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═e2f50ace-ad14-4e7c-af74-abf3ca9df9fb
# ╠═35da0c45-6d91-4ad9-811f-5d633684208e
# ╠═c3a5e97e-00ed-44af-9c67-fa8198900fbd
# ╠═82c0d6f9-373d-4866-81eb-9cd2d0981310
# ╠═2e6bcf9b-002f-4fbc-80c1-0cfd2ab1d253
# ╠═4c06d21c-ac14-4522-bf25-2e0a1ed2d6b9
# ╠═1e3da0c9-2f96-4637-9093-ac7f10c1ad27
# ╠═48c73104-fe4c-4543-942a-6f23b0fd2547
# ╠═48c73104-fe4c-4543-0001-6f23b0fd2547
# ╠═48c73104-fe4c-4543-0002-6f23b0fd2547
# ╠═c4402f72-67ac-4630-a651-da81c1df71bf
# ╠═e796d5be-e011-45d3-ad90-58769feb5e85
# ╠═d890aff1-dbd0-451c-bf14-bde9758c3be0
# ╠═fdbbc8d9-83d6-413e-aff8-aca18f24dfea
# ╠═cc0cae75-ba10-4a62-b0ef-22259e40a083
# ╠═855553e3-491b-4e8a-a482-95855697e063
# ╠═1f7a7093-ce8e-461f-8b91-69266de86748
# ╠═f83fc65f-5f7b-498d-9ed3-0762565ad710
# ╠═4714c6ae-27d9-47db-b12e-126283b10606
# ╠═f9cc2880-bab1-4be3-b67b-d43508df8d3b
# ╠═1f2b9b74-de46-401d-8a46-0434b9f9aca1
# ╠═6a5cb372-60cb-4ffc-b4f0-22e4016104e7
# ╠═e832c1b7-8c04-4146-90f6-1628e91fea2a
# ╠═64cf1b8b-6686-4eeb-a3cc-786300ea7c7d
# ╠═1d023a0c-a93a-451c-a894-1d1f6a4b78a9
# ╠═b8ce87d4-4768-4e8a-a16e-9e68b00b6617
# ╠═e3c221c6-6c4a-4b5f-93a6-30d3508ac9d2
# ╠═8368ae01-ce53-449e-87dd-8dfa3f29f8f4
# ╠═f207d849-ebff-4e6c-95bb-50693cb7c9b6
# ╠═3033d5c1-d8e0-4d46-0001-7dec4ff7afbd
# ╠═3033d5c1-d8e0-4d46-a4b9-7dec4ff7afbd
# ╠═528beffe-0707-4661-8970-def1b1e00ea5
# ╠═f951805e-515a-475f-893f-bb8b968e425c
# ╠═86ae20a9-e69c-4d63-9119-395449e9ac09
# ╠═9be2fd5c-4b6c-4e13-b2aa-fb7120a504b7
# ╠═16782c72-ecb1-48ec-8510-78e2e0689a10
# ╠═7344190c-7989-4b55-b7be-357f7d6b7370
# ╠═11938cb6-46b3-0001-96c0-ef6424d1d0db
# ╟─11938cb6-46b3-0002-96c0-ef6424d1d0db
# ╠═11938cb6-46b3-0003-96c0-ef6424d1d0db
# ╠═11938cb6-46b3-499b-96c0-ef6424d1d0db
# ╠═d2c248b1-c48e-437b-a910-edcc59b4424f
# ╠═7b306ed5-4bda-465d-abf2-4d07cb4642c1
# ╠═c45f6aac-bff3-4c98-8bd5-c91a98c9eef7
# ╠═1f1505a6-c6af-4fdd-9583-6e783e76de4f
# ╠═980e48bd-9dd1-4195-8086-40785a7f43e1
# ╠═fc146750-46b2-4084-a6d2-0d91c0e104e6
# ╠═fc3e96a6-7af8-4ef0-a4ca-b33509aa512f
# ╠═ed2f2fa3-be19-4e2d-ad28-36cb053ed6bd
# ╠═6151e4d1-a7a6-4fa7-bc77-0e7f3b2a3cc0
# ╠═b5449c6a-9c43-4f9c-81f7-f01277acb109
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
