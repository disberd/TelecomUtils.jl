### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 13646410-4e96-11ec-3e3d-99763ba1aeea
begin
	using CoordinateTransformations
	using StaticArrays
	using LinearAlgebra
	using Unitful
	using Unitful.DefaultSymbols
	using Rotations
	using Parameters
	using SatelliteToolbox: geodetic_to_ecef, ecef_to_geodetic, wgs84_ellipsoid
	using DocStringExtensions
	using Parameters
	import Proj4
end

# ╔═╡ 73e00cef-9734-439a-b89b-7c1d99aab74e
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	using BenchmarkTools
	using PlutoTest
	using PlutoUtils
	using PlotlyBase
	using PlutoDevMacros
end
  ╠═╡ =#

# ╔═╡ 3bda5426-c0de-493f-9514-30b6fe762463
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Packages
"""
  ╠═╡ =#

# ╔═╡ dcc81988-903b-4707-a70c-09c38682c80f
# ╠═╡ skip_as_script = true
#=╠═╡
ToC()
  ╠═╡ =#

# ╔═╡ 3fd1046c-fabf-4264-9638-ba41301b1804
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## load other notebook
"""
  ╠═╡ =#

# ╔═╡ 7729ce27-df74-4393-ab70-c4e2864c85f5
@plutoinclude "satview_transformations.jl" "all"

# ╔═╡ de735c56-612c-4ffd-8335-95f20a129390
# ╠═╡ skip_as_script = true
#=╠═╡
@macroexpand @plutoinclude "satview_transformations.jl" "all"
  ╠═╡ =#

# ╔═╡ a57e3983-21de-4a2e-a227-8265fee6b56b
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Exports
"""
  ╠═╡ =#

# ╔═╡ 030e15c5-83a8-4a24-836a-96b6f4f0bb04
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# SatView
"""
  ╠═╡ =#

# ╔═╡ b966fa0c-dc52-4821-bc32-e78dd3272ce1
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
We want to define a SatView mutable struct that represents a satellite and can be used to compute view angles, ranges and ERA for specific points on ground
"""
  ╠═╡ =#

# ╔═╡ e469d967-79e4-4ef2-b635-51a183cb12e7
begin
	"""
	SatView(lla::LLA,earthmodel::EarthModel)
	SatView(ecef::StaticVector{3},earthmodel::EarthModel)
Object representing the instantaneous position of a satellite and used to compute various view angles to and from points on ground, as well as inverse geodesic computations.

# Fields
$TYPEDFIELDS

If multiple satellites have to be tracked, the `EarthModel` instance `earthmodel` should be generated once and then passed to all the SatView instances to ensure that all the satellites are referring to the same earth model.\\
Doing so will enforce the same Ellipsoid is shared between all satellites even when it is changed from one of the SatView instances.

See also: [`change_position!`](@ref), [`get_range`](@ref), [`get_pointing`](@ref), [`get_lla`](@ref), [`get_ecef`](@ref), [`geod_inverse`](@ref), [`get_distance_on_earth`](@ref), [`get_nadir_beam_diameter`](@ref), [`ExtraOutput`](@ref).
	"""
	mutable struct SatView
		"ECEF coordinates of the current satellite position"
		ecef::SVector{3,Float64}
		"LLA coordinates of the current satellite position"
		lla::LLA
		"Reference EarthModel used for the projections and for the geodesic computations"
		earthmodel::EarthModel
		"Rotation Matrix to go from the nadir-pointing CRS to the ECEF CRS"
		R::RotMatrix3{Float64}
	end

	# Custom constructor
	function SatView(lla::LLA,em::EarthModel)
		ecef = ECEFfromLLA(em.ellipsoid)(lla)
		R = _rotation_matrix(:ECEFfromUV, lla.lat, lla.lon)
		SatView(ecef,lla,em,R)
	end
	function SatView(ecef::StaticVector{3},em::EarthModel)
		lla = LLAfromECEF(em.ellipsoid)(ecef)
		R = _rotation_matrix(:ECEFfromUV, lla.lat, lla.lon)
		SatView(ecef,lla,em,R)
	end
end

# ╔═╡ 4794b9dc-4297-402b-94a2-ed686584bb09
function Base.getproperty(sv::SatView, name::Symbol)
	if name ∈ (:ellipsoid, :geod)
		return getfield(getfield(sv,:earthmodel),name)
	else
		return getfield(sv, name)
	end
end

# ╔═╡ dd8d119d-550c-4217-923d-43aaf2b8327b
function Base.setproperty!(sv::SatView, name::Symbol, x)
	if name ∈ (:ellipsoid, :geod)
		return setproperty!(getfield(sv,:earthmodel),name, x)
	else
		return setfield!(sv, name, x)
	end
end

# ╔═╡ 6eb9424e-3dd2-46d4-b4d2-81596bb81668
# ╠═╡ skip_as_script = true
#=╠═╡
em = EarthModel()
  ╠═╡ =#

# ╔═╡ 170f2914-fdf9-46c8-a8e0-9130b046bd60
# ╠═╡ skip_as_script = true
#=╠═╡
sv = SatView(LLA(0,0,100km), em)
  ╠═╡ =#

# ╔═╡ e28cbfc7-408e-49b5-9aeb-bd01c32fba46
# ╠═╡ skip_as_script = true
#=╠═╡
let
	sv = SatView(LLA(0,0,100km), EarthModel())
	sv.ellipsoid = wgs84_ellipsoid
	sv.geod
end
  ╠═╡ =#

# ╔═╡ 5f1fd82d-f441-4a1b-9840-773a8635d3db
# ╠═╡ skip_as_script = true
#=╠═╡
@benchmark getproperty($sv, :geod)
  ╠═╡ =#

# ╔═╡ 091e4ec2-ea9e-411e-8f39-73aeb73c0214
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Change Position
"""
  ╠═╡ =#

# ╔═╡ 41370c82-d32a-41ea-a21a-614574292c21
begin
	"""
	change_position!(sv::SatView, ecef, lla::LLA, R)
	change_position!(sv::SatView, lla::LLA)
	change_position!(sv::SatView, ecef::StaticVector{3})
Change the position of a [`SatView`](@ref) object `sv`, also returning as output the modified
`sv`.  The function mutates the `ecef`, `lla` and `R` fields of the `sv`
object with the values provided as arguments (when using the 1st method above).\\
If only `ecef` or `lla` coordinates are provided (2nd and 3rd method above), the remaining
two arguments are computed automatically.

One would normally use either the 2nd or 3rd method so that the two missing components are
correctly computed by the function.\\
The first method avoids computations but does not validate that `ecef`, `lla` and `R` are
correct and refer to the same position in space. For this reason the first method should
only be used if those values are correctly pre-computed elsewhere and one wants to avoid the
duplicate computations. 

See also: [`SatView`](@ref), [`get_range`](@ref), [`get_pointing`](@ref), [`get_lla`](@ref),
[`get_ecef`](@ref), [`get_ecef`](@ref).
	"""
	function change_position!(sv::SatView, ecef, lla::LLA, R)
		setfield!(sv,:ecef,ecef)
		setfield!(sv,:lla,lla)
		setfield!(sv,:R,R)
		return sv
	end
	function change_position!(sv::SatView, lla::LLA)
		ecef = ECEFfromLLA(sv.ellipsoid)(lla)
		R = _rotation_matrix(:ECEFfromUV, lla.lat, lla.lon)
		change_position!(sv, ecef, lla, R)
	end
	function change_position!(sv::SatView, ecef::StaticVector{3})
		lla = LLAfromECEF(sv.ellipsoid)(ecef)
		R = _rotation_matrix(:ECEFfromUV, lla.lat, lla.lon)
		change_position!(sv, ecef, lla, R)
	end
end

# ╔═╡ f77dcbcd-f042-4f7c-b97f-de63637229d0
# ╠═╡ skip_as_script = true
#=╠═╡
@benchmark change_position!($sv,$(SA_F64[1e7,0,0]))
  ╠═╡ =#

# ╔═╡ 78a8e7a4-333d-44ca-a438-fd85d7078300
# ╠═╡ skip_as_script = true
#=╠═╡
@code_warntype change_position!(sv,SA_F64[1e7,0,0])
  ╠═╡ =#

# ╔═╡ 709c44c8-c580-4cf7-8376-9c513eb3bd53
# ╠═╡ skip_as_script = true
#=╠═╡
@benchmark change_position!($sv,LLA(0,0,100km))
  ╠═╡ =#

# ╔═╡ 84769564-8ba8-46f5-b494-b0689d9abd65
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Get Range
"""
  ╠═╡ =#

# ╔═╡ 642f2ede-b154-4260-a959-0a47ca4793b7
"""
	get_range(sv::SatView,uv::StaticVector{2,<:Number},h::Real = 0.0)
	get_range(sv::SatView,uv::Tuple{<:Number,<:Number},h::Real = 0.0)
Get the range [in m] between the satellite and a given point identified by the uv pointing and the reference altitude above the ellipsoid.
"""
function get_range(sv::SatView,uv::StaticVector{2,<:Number},h::Real = 0.0)
	# Find the ecef coordinate of the target point on earth
	ecef = ECEFfromUV(sv.ecef, sv.R, sv.ellipsoid)(uv,h)
	# Return the distance betwteen the satellite and the point
	return norm(ecef-sv.ecef)
end

# ╔═╡ 7a75d351-8583-455f-89c4-2d50cf79ea96
get_range(sv::SatView,uv::Tuple{<:Number,<:Number},args...) = get_range(sv,SVector(uv),args...)

# ╔═╡ 556934d8-d3ee-4a43-8f74-0939c5431c6f
"""
	get_range(sv::SatView,lla::LLA)
	get_range(sv::SatView,lla::LLA, ::ExtraOutput)
	get_range(sv::SatView,ecef::StaticVector{3})
	get_range(sv::SatView,ecef::StaticVector{3}, ::ExtraOutput)
Get the range [in m] between the satellite and a given point identified either by its LLA or ECEF coordinates. Returns NaN if the point is not visible either because it is covered by the earth surface or it is *above* the satellite.

If called with an instance of the `ExtraOutput` struct as last argument, it also provides the WND coordinates of the target point as seen from the satellite

See also: [`SatView`](@ref), [`change_position!`](@ref), [`get_pointing`](@ref), [`get_lla`](@ref), [`get_ecef`](@ref).
"""
function get_range(sv::SatView, ecef::StaticVector{3}, ::ExtraOutput)
	Δecef = ecef - sv.ecef
	dist = norm(Δecef)
	# Find if the target point is below the satellite, we do this by checking the last coordinate of the WND coordinates of the point
	wnd = sv.R'Δecef
	wnd[3] < 0 && return NaN, wnd
	# We have to check that the given lla is visible from the satellite, this happens if either there is no intersection with earth in the direction of pointing, or if the first intersection happens for a range greater than th computed one
	t₁, t₂ = _intersection_solutions(Δecef/dist, sv.ecef, sv.ellipsoid.a, sv.ellipsoid.b)
	r = (isnan(t₁) || t₁ > dist-1e-5) ? dist : NaN
	return r,wnd
end

# ╔═╡ 68417643-fa77-4780-9890-b0dac95bdb7f
get_range(sv::SatView, ecef::StaticVector{3}) = get_range(sv,ecef,ExtraOutput())[1]

# ╔═╡ da78f52b-30b6-4faf-bcea-b665c10ff4fe
get_range(sv::SatView, lla::LLA, args...) = get_range(sv,ECEFfromLLA(sv.ellipsoid)(lla),args...)

# ╔═╡ 449b49de-2951-41fc-ba46-89eaa6c52e79
# ╠═╡ skip_as_script = true
#=╠═╡
get_range(sv,LLA(0°,5°,10km),ExtraOutput())
  ╠═╡ =#

# ╔═╡ e7443f5b-a1a8-4866-9a64-ce7587465911
# ╠═╡ skip_as_script = true
#=╠═╡
get_range(SatView(LLA(0,0,600km),em),(0,0))
  ╠═╡ =#

# ╔═╡ 7c07a3c1-c1ec-4b83-b7c6-251edf91273c
# ╠═╡ skip_as_script = true
#=╠═╡
@benchmark $get_range($sv,LLA(0,0,10km))
  ╠═╡ =#

# ╔═╡ 39a1850b-f64a-4157-8f07-d7a78918fea1
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Get Pointing
"""
  ╠═╡ =#

# ╔═╡ 51987c04-18f5-46bb-a3ba-5f94907a7960
"""
	get_pointing(sv::SatView,lla::LLA,kind::Symbol=:uv)
	get_pointing(sv::SatView,ecef::StaticVector{3},kind::Symbol=:uv)
Provide the 2-D angular pointing at which the target point (specified as LLA or ECEF) is seen from the satellite.

`kind` is used to select whether the output should be given in UV or ThetaPhi coordinates. The result is provided as ThetaPhi [in rad] if `kind ∈ (:ThetaPhi, :thetaphi, :θφ)`

See also: [`SatView`](@ref), [`get_range`](@ref), [`get_pointing`](@ref), [`get_lla`](@ref), [`get_ecef`](@ref), [`get_distance_on_earth`](@ref).
"""
function get_pointing(sv::SatView,ecef::StaticVector{3},kind::Symbol=:uv)
	uv = UVfromECEF(sv.ecef,sv.R',sv.ellipsoid)(ecef)
	if kind ∈ (:ThetaPhi, :thetaphi, :θφ)
		return ThetaPhifromUV()(uv)
	else
		return uv
	end
end

# ╔═╡ a6db34bc-b846-49aa-8d57-fb32cdce1684
get_pointing(sv::SatView,lla::LLA,args...) = get_pointing(sv,ECEFfromLLA(sv.ellipsoid)(lla),args...)

# ╔═╡ 1758748c-fa4b-4414-a05d-a32970c7a94b
# ╠═╡ skip_as_script = true
#=╠═╡
get_pointing(SatView(LLA(0,0,600km),em),LLA(1°,1°,0),:thetaphi)
  ╠═╡ =#

# ╔═╡ cc1c1137-a253-49de-8293-5819236a00cf
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Get LLA
"""
  ╠═╡ =#

# ╔═╡ 1f7bf45c-b33b-4bfe-b82d-05b908ce375e
"""
	get_lla(sv::SatView,pointing::Point2D,kind::Symbol=:uv; h = 0.0)
Computes the LLA coordinates of a point on earth seen from the satellite identified by `sv`, the angular pointing identified by `pointing` and located at an altitude `h` [m] above the reference Earth ellipsoid of `sv`.

The optional argument `kind` is used to select whether the pointing is expressed in ThetaPhi (`kind ∈ (:ThetaPhi, :thetaphi, :θφ)`) [rad] or UV coordinates.

See also: [`SatView`](@ref), [`get_range`](@ref), [`get_pointing`](@ref), [`get_lla`](@ref), [`get_ecef`](@ref), [`get_distance_on_earth`](@ref).
"""
function get_lla(sv::SatView,pointing::Point2D,kind::Symbol=:uv; h = 0.0)
	uv = if kind ∈ (:ThetaPhi, :thetaphi, :θφ)
		UVfromThetaPhi()(pointing)
	else
		pointing
	end
	lla = LLAfromUV(sv.ecef,sv.R,sv.ellipsoid)(uv)
end

# ╔═╡ 1f27b72f-9a3b-4732-a98e-d216af067072
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Get ECEF
"""
  ╠═╡ =#

# ╔═╡ 948cc7a1-d85e-4cfe-b2e4-e047bcbac305
"""
	get_ecef(sv::SatView,pointing::Point2D,kind::Symbol=:uv; h = 0.0)
Computes the ECEF coordinates of a point on earth seen from the satellite identified by `sv`, the angular pointing identified by `pointing` and located at an altitude `h` [m] above the reference Earth ellipsoid of `sv`.

The optional argument `kind` is used to select whether the pointing is expressed in ThetaPhi (`kind ∈ (:ThetaPhi, :thetaphi, :θφ)`) [rad] or UV coordinates.

See also: [`SatView`](@ref), [`get_range`](@ref), [`get_pointing`](@ref), [`get_lla`](@ref), [`change_position!`](@ref), [`get_distance_on_earth`](@ref).
"""
function get_ecef(sv::SatView,pointing::Point2D,kind::Symbol=:uv; h = 0.0)
	uv = if kind ∈ (:ThetaPhi, :thetaphi, :θφ)
		UVfromThetaPhi()(pointing)
	else
		pointing
	end
	lla = ECEFfromUV(sv.ecef,sv.R,sv.ellipsoid)(uv)
end

# ╔═╡ 8bc60d8d-7b54-4dce-a3e4-e336c0b16d4e
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Get ERA
"""
  ╠═╡ =#

# ╔═╡ b2cb0afd-1220-40bd-8e1b-6df35e3db2f1
function get_era(sv::SatView, lla::LLA)
	ERAfromECEF(lla;ellipsoid = sv.ellipsoid)(sv.ecef)
end

# ╔═╡ 8fccb117-2048-4607-8db1-f8df7f5ef156
function get_era(sv::SatView, ecef::StaticVector{3})
	ERAfromECEF(ecef;ellipsoid = sv.ellipsoid)(sv.ecef)
end

# ╔═╡ ee657a11-c976-4128-8bb4-2336a5ecd319
# ╠═╡ skip_as_script = true
#=╠═╡
# We test that a non-visible point is nan
@test get_era(SatView(LLA(0,0,700km),em),LLA(40°,-39°,0)) |> isnan
  ╠═╡ =#

# ╔═╡ 2ad13505-0c60-4ccb-b536-e865c24a0396
# ╠═╡ skip_as_script = true
#=╠═╡
# We test that a visible point is not nan
@test get_era(SatView(LLA(0,0,600km),em),LLA(0,0,500km)) ≈ ERA(90°, 100km, 0°)
  ╠═╡ =#

# ╔═╡ 97c3ab73-5d2b-4871-aaa2-f8d7f1a7204d
# ╠═╡ skip_as_script = true
#=╠═╡
@benchmark get_era($sv,LLA(0,0,0))
  ╠═╡ =#

# ╔═╡ 80b5256d-e3f1-4329-be97-34e557377466
function meshgrid(xin,yin)
	nx=length(xin)
	ny=length(yin)
	xout=zeros(ny,nx)
	yout=zeros(ny,nx)
	for jx=1:nx
	    for ix=1:ny
	        xout[ix,jx]=xin[jx]
	        yout[ix,jx]=yin[ix]
	    end
	end
	return xout,yout
end;

# ╔═╡ bbf6f990-40b3-471f-a46c-61f5fd6f5824
# Visual Test ERA
begin
	earthmodel = EarthModel()
	svTest = SatView(LLA(0,0,700km),earthmodel)
	gridRes = 0.5
	x_plot = -180:gridRes:180
	y_plot = -90:gridRes:90

	lonMesh,latMesh = meshgrid(-180:gridRes:180,-90:gridRes:90)
	el_plot = fill(NaN,size(lonMesh,1),size(lonMesh,2))
	az_plot = fill(NaN,size(lonMesh,1),size(lonMesh,2))
	r_plot = fill(NaN,size(lonMesh,1),size(lonMesh,2))
	for row = 1:size(latMesh,1)
		for col = 1:size(latMesh,2)
			era = get_era(svTest,LLA(latMesh[row,col]*π/180,lonMesh[row,col]*π/180,0))
			el_plot[row,col] = rad2deg(era.el)
			az_plot[row,col] = rad2deg(era.az)
			r_plot[row,col] = era.r
		end
	end
end

# ╔═╡ 74b99a07-9e73-4532-a50d-10221c47f324
# ╠═╡ skip_as_script = true
#=╠═╡
let
	el_plot = heatmap(
		x = x_plot,
		y = y_plot,
		z = el_plot
	)
	
	Plot(el_plot, Layout(
					yaxis_title = "LAT",
					xaxis_title = "LON",
					title = "Elevation Test"))
end
  ╠═╡ =#

# ╔═╡ 1c0aa81c-efa2-4ba4-a3b2-70276d76c4f1
# ╠═╡ skip_as_script = true
#=╠═╡
let
	az_plot = heatmap(
		x = x_plot,
		y = y_plot,
		z = az_plot
	)
	
	Plot(az_plot, Layout(
					yaxis_title = "LAT",
					xaxis_title = "LON",
					title = "Azimuth Test"))
end
  ╠═╡ =#

# ╔═╡ 5023b71d-219e-4f2f-b319-e9899e9702ac
# ╠═╡ skip_as_script = true
#=╠═╡
let
	r_plot = heatmap(
		x = x_plot,
		y = y_plot,
		z = r_plot
	)
	
	Plot(r_plot, Layout(
					yaxis_title = "LAT",
					xaxis_title = "LON",
					title = "Range Test"))
end
  ╠═╡ =#

# ╔═╡ 64370881-a469-4748-97c5-ec27199d529b
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Get Distance on Earth
"""
  ╠═╡ =#

# ╔═╡ 2efb01b8-16b1-4186-94f4-cdfbca1310de
geod_inverse(sv::SatView, args...) = geod_inverse(sv.geod,args...)

# ╔═╡ e0915eab-a53d-4fb2-9029-83793073ac3c
"""
	get_distance_on_earth(sv::SatView, p1::Point2D, p2::Point2D, kind::Symbol=:uv)
	get_distance_on_earth(sv::SatView, p1::Point2D, p2::Point2D, ::ExtraOutput, kind::Symbol=:uv)
	get_distance_on_earth(sv::SatView, lla1::LLA, lla2::LLA)
	get_distance_on_earth(sv::SatView, lla1::LLA, lla2::LLA, ::ExtraOutput)
Computes the distance [m] between the points on the earth surface (`lla1` and `lla2`) using the reference earth model used by `sv`.

If the points are not provided as LLA instances, but as angular directions (`p1` and `p2`), `lla1` and `lla2` as first computed from `p1` and `p2` using the SatView object `sv` as reference.

When called with angular directions, the optional argument `kind` is used to select whether the pointing is expressed in ThetaPhi (`kind ∈ (:ThetaPhi, :thetaphi, :θφ)`) [rad] or UV coordinates.

If an instance of `ExtraOutput` is provided as 4th argument, the function also returns the (forward) azimuth angle between `lla1` and `lla2` (2nd output, [deg]) and the azimuth angle between `lla2` and `lla1` (third output, [deg])

See also: [`SatView`](@ref), [`get_range`](@ref), [`get_pointing`](@ref), [`get_lla`](@ref), [`get_ecef`](@ref), [`geod_inverse`](@ref), [`get_nadir_beam_diameter`](@ref), [`ExtraOutput`](@ref).
"""
function get_distance_on_earth(sv::SatView, p1::Point2D, p2::Point2D, eo::ExtraOutput, kind::Symbol=:uv)
	lla1 = get_lla(sv, p1, kind)
	lla2 = get_lla(sv, p2, kind)
	get_distance_on_earth(sv, lla1, lla2, eo)
end

# ╔═╡ 30d32b7a-c95c-4d80-a60a-a87b27b3bf3c
get_distance_on_earth(sv::SatView, p1::Point2D, p2::Point2D, kind::Symbol=:uv) = get_distance_on_earth(sv, p1, p2, ExtraOutput(), kind)[1]

# ╔═╡ 407101b2-c794-49b0-9f8b-07fb45b80ca9
get_distance_on_earth(sv::SatView, lla1::LLA, lla2::LLA, ::ExtraOutput) = geod_inverse(sv.geod, lla1, lla2)

# ╔═╡ af71267d-b5ee-46b7-bf8d-d740033d35e0
get_distance_on_earth(sv::SatView, lla1::LLA, lla2::LLA) = get_distance_on_earth(sv, lla1, lla2, ExtraOutput())[1]

# ╔═╡ 2612961b-0e1e-4615-8959-74ab3bc919f9
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Get Nadir Beam Diameter
"""
  ╠═╡ =#

# ╔═╡ 30959832-9eb2-48c5-83d5-776d336c9aa7
"""
$SIGNATURES
Computes the diameter [m] on earth of a beam pointed at nadir from the satellite position identified by `sv` and assuming a 3db beamwidth identified by the scan angle `scan_3db` [deg] at which the beam pattern is 3dB below the peak.

The computation computes the diameter along the `U` direction and takes into account the reference ellipsoid of `sv`, so the resulting diameter is dependent on the satellite lat/lon position

See also: [`SatView`](@ref), [`get_range`](@ref), [`get_pointing`](@ref), [`get_lla`](@ref), [`get_ecef`](@ref), [`geod_inverse`](@ref), [`get_distance_on_earth`](@ref). 
"""
function get_nadir_beam_diameter(sv, scan_3db)
	uv_radius = sind(scan_3db)
	p1 = (uv_radius, 0)
	p2 = (-uv_radius, 0)
	get_distance_on_earth(sv, p1, p2)
end

# ╔═╡ 4af7a092-8f42-4aef-9c09-feab8ebc1d87
# ╠═╡ skip_as_script = true
#=╠═╡
get_nadir_beam_diameter(SatView(LLA(50°,0°,735km), EarthModel()), 55)
  ╠═╡ =#

# ╔═╡ 22a1312d-c7f4-486f-a62c-3389f7715bae
md"""
## Get Angle Between Satellites
"""

# ╔═╡ 28938135-9c93-4992-bfec-f106d8aa0bf6
function get_angle2Sat(sv1::SatView, sv2::SatView, ::ExtraOutput)  
	# Find the angle between 2 satellites (same or different shell) from the nadir/zenith of sv1 to sv2 and viceversa. If the LoS between the two satellites is obstructed by Earth return NaN.

	# Find the difference vector between the satellites (vector "connecting" sv1 and sv2)
 	pdiff = (sv1.ecef - sv2.ecef) 
 	 
 	# Find the magnitude of the difference to compare with the intersection solutions 
 	t = norm(pdiff) 
 	 
 	# Find the intersection points with the ellipsoid 
 	t₁,t₂ = _intersection_solutions(pdiff./t, sv1.ecef, sv1.earthmodel.ellipsoid.a, sv1.earthmodel.ellipsoid.b) 
 	 
 	# If both t₁ and t₂ are NaN, it means that no intersection with the ellipsoid is found and so there is no earth blockage 
 	# If t <= t₁ also no blockage is present 
 	# If t > t₁ then the earth is blocking the view point so we return NaN 
 	 
 	# The 1e-3 is there because the computed distance might have some error that is usually way below one mm, and 1mm shouldn't change anything for our required precision 
	# Return NaN if the Earth is blocking the LoS between 2 satellites
 	!isnan(t₁) && abs(t) > abs(t₁)+1e-3 && return NaN, NaN, NaN 
 	 
 	# Find the coordinates in the West-North-Down CRS (centered in sv1)
 	wnd = sv1.R * pdiff 

	# Convert in spherical coordinates
	x,y,z = wnd
 	r = hypot(x, y, z) 
 	θ = r == 0 ? 0 : acos(z/r) # 0: nadir
	# //TODO: check consistency with the rest of the cde for angle conventions (remove -x)
 	ϕ = r == 0 ? 0 : atan(y,x) # angle measured from West to North clockwise wrt the local reference 
 	 
 	# Return coordinates 
	# //TODO: check unit for consistency
 	# return θ * rad, ϕ * rad, r * m
 	return θ, ϕ, r
 end

# ╔═╡ 56d88bb9-3b33-4b1a-88ae-d90af4de2bd1
# Call returning only angle from sv1 to sv2
get_angle2Sat(sv1::SatView, sv2::SatView) = get_angle2Sat(sv1,sv2,ExtraOutput())[1]

# ╔═╡ b9dacaaf-b55c-46c8-8fd0-ad520505ecbb
export SatView, change_position!, get_range, get_era, get_pointing, get_lla, get_ecef, get_distance_on_earth, get_nadir_beam_diameter, get_angle2Sat

# ╔═╡ c02d0705-6647-4a44-8ae8-fc256f18c4ce
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Tests
"""
  ╠═╡ =#

# ╔═╡ d15726ab-5a28-4a24-b5ed-b3c8ecb6c581
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## nadir beam diameter
"""
  ╠═╡ =#

# ╔═╡ ae686da9-45d5-4fc2-9cbd-2d828d792407
# ╠═╡ skip_as_script = true
#=╠═╡
@test get_nadir_beam_diameter(SatView(LLA(90°,0°,735km), EarthModel()), 55) ≈ get_nadir_beam_diameter(SatView(LLA(0°,0°,735km), EarthModel()), 55)
  ╠═╡ =#

# ╔═╡ 71d3f92e-d143-40dc-8701-37f9053766ef
# ╠═╡ skip_as_script = true
#=╠═╡
@test get_nadir_beam_diameter(SatView(LLA(90°,0°,735km), EarthModel(wgs84_ellipsoid)), 55) ≉ get_nadir_beam_diameter(SatView(LLA(0°,0°,735km), EarthModel(wgs84_ellipsoid)), 55)
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
CoordinateTransformations = "150eb455-5306-5404-9cee-2592286d6298"
DocStringExtensions = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
PlotlyBase = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
PlutoDevMacros = "a0499f29-c39b-4c5c-807c-88074221b949"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUtils = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
Proj4 = "9a7e659c-8ee8-5706-894e-f68f43bc57ea"
Rotations = "6038ab10-8711-5258-84ad-4b1120ba62dc"
SatelliteToolbox = "6ac157d9-b43d-51bb-8fab-48bf53814f4a"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
BenchmarkTools = "~1.3.1"
CoordinateTransformations = "~0.6.2"
DocStringExtensions = "~0.8.6"
Parameters = "~0.12.3"
PlotlyBase = "~0.8.18"
PlutoDevMacros = "~0.4.5"
PlutoTest = "~0.2.2"
PlutoUtils = "~0.5.9"
Proj4 = "~0.7.6"
Rotations = "~1.2.0"
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

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "4c10eee4af024676200bc7752e536f858c6b8f93"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.1"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.Chain]]
git-tree-sha1 = "8c4920235f6c561e401dfe569beb8b924adad003"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.5.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9489214b993cd42d17f44c36e359bf6a7c919abf"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "1e315e3f4b0b7ce40feded39c73049692126cf53"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.3"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

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

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

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

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "b3364212fb5d870f724876ffcd34dd8ec6d98918"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.7"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

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

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

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

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "09e4b894ce6a976c354a69041a04748180d43637"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.15"

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

[[deps.NaNMath]]
git-tree-sha1 = "737a5957f387b17e74d4ad2f440eb330b39a62c5"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.0"

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

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

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

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "180d744848ba316a3d0fdf4dbd34b77c7242963a"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.18"

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

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.Proj4]]
deps = ["CEnum", "CoordinateTransformations", "PROJ_jll", "StaticArrays"]
git-tree-sha1 = "5f15f1c647b563e49f655fbbfd4e2ade24bd3c64"
uuid = "9a7e659c-8ee8-5706-894e-f68f43bc57ea"
version = "0.7.6"

[[deps.Quaternions]]
deps = ["DualNumbers", "LinearAlgebra", "Random"]
git-tree-sha1 = "0b345302b17b0e694092621915de0e0dc7443a1a"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.4.9"

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

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "405148000e80f70b31e7732ea93288aecb1793fa"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.2.0"

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

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "a9e798cae4867e3a41cae2dd9eb60c047f1212db"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.6"

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

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

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
# ╠═3bda5426-c0de-493f-9514-30b6fe762463
# ╠═13646410-4e96-11ec-3e3d-99763ba1aeea
# ╠═73e00cef-9734-439a-b89b-7c1d99aab74e
# ╠═dcc81988-903b-4707-a70c-09c38682c80f
# ╠═3fd1046c-fabf-4264-9638-ba41301b1804
# ╠═7729ce27-df74-4393-ab70-c4e2864c85f5
# ╠═de735c56-612c-4ffd-8335-95f20a129390
# ╠═a57e3983-21de-4a2e-a227-8265fee6b56b
# ╠═b9dacaaf-b55c-46c8-8fd0-ad520505ecbb
# ╠═030e15c5-83a8-4a24-836a-96b6f4f0bb04
# ╠═b966fa0c-dc52-4821-bc32-e78dd3272ce1
# ╠═e469d967-79e4-4ef2-b635-51a183cb12e7
# ╠═4794b9dc-4297-402b-94a2-ed686584bb09
# ╠═dd8d119d-550c-4217-923d-43aaf2b8327b
# ╠═6eb9424e-3dd2-46d4-b4d2-81596bb81668
# ╠═170f2914-fdf9-46c8-a8e0-9130b046bd60
# ╠═e28cbfc7-408e-49b5-9aeb-bd01c32fba46
# ╠═5f1fd82d-f441-4a1b-9840-773a8635d3db
# ╠═091e4ec2-ea9e-411e-8f39-73aeb73c0214
# ╠═41370c82-d32a-41ea-a21a-614574292c21
# ╠═f77dcbcd-f042-4f7c-b97f-de63637229d0
# ╠═78a8e7a4-333d-44ca-a438-fd85d7078300
# ╠═709c44c8-c580-4cf7-8376-9c513eb3bd53
# ╠═84769564-8ba8-46f5-b494-b0689d9abd65
# ╠═642f2ede-b154-4260-a959-0a47ca4793b7
# ╠═7a75d351-8583-455f-89c4-2d50cf79ea96
# ╠═556934d8-d3ee-4a43-8f74-0939c5431c6f
# ╠═68417643-fa77-4780-9890-b0dac95bdb7f
# ╠═da78f52b-30b6-4faf-bcea-b665c10ff4fe
# ╠═449b49de-2951-41fc-ba46-89eaa6c52e79
# ╠═e7443f5b-a1a8-4866-9a64-ce7587465911
# ╠═7c07a3c1-c1ec-4b83-b7c6-251edf91273c
# ╠═39a1850b-f64a-4157-8f07-d7a78918fea1
# ╠═51987c04-18f5-46bb-a3ba-5f94907a7960
# ╠═a6db34bc-b846-49aa-8d57-fb32cdce1684
# ╠═1758748c-fa4b-4414-a05d-a32970c7a94b
# ╠═cc1c1137-a253-49de-8293-5819236a00cf
# ╠═1f7bf45c-b33b-4bfe-b82d-05b908ce375e
# ╠═1f27b72f-9a3b-4732-a98e-d216af067072
# ╠═948cc7a1-d85e-4cfe-b2e4-e047bcbac305
# ╠═8bc60d8d-7b54-4dce-a3e4-e336c0b16d4e
# ╠═b2cb0afd-1220-40bd-8e1b-6df35e3db2f1
# ╠═8fccb117-2048-4607-8db1-f8df7f5ef156
# ╠═ee657a11-c976-4128-8bb4-2336a5ecd319
# ╠═2ad13505-0c60-4ccb-b536-e865c24a0396
# ╠═97c3ab73-5d2b-4871-aaa2-f8d7f1a7204d
# ╠═bbf6f990-40b3-471f-a46c-61f5fd6f5824
# ╠═74b99a07-9e73-4532-a50d-10221c47f324
# ╠═1c0aa81c-efa2-4ba4-a3b2-70276d76c4f1
# ╠═5023b71d-219e-4f2f-b319-e9899e9702ac
# ╠═80b5256d-e3f1-4329-be97-34e557377466
# ╠═64370881-a469-4748-97c5-ec27199d529b
# ╠═2efb01b8-16b1-4186-94f4-cdfbca1310de
# ╠═e0915eab-a53d-4fb2-9029-83793073ac3c
# ╠═30d32b7a-c95c-4d80-a60a-a87b27b3bf3c
# ╠═407101b2-c794-49b0-9f8b-07fb45b80ca9
# ╠═af71267d-b5ee-46b7-bf8d-d740033d35e0
# ╟─2612961b-0e1e-4615-8959-74ab3bc919f9
# ╠═30959832-9eb2-48c5-83d5-776d336c9aa7
# ╠═4af7a092-8f42-4aef-9c09-feab8ebc1d87
# ╟─22a1312d-c7f4-486f-a62c-3389f7715bae
# ╠═28938135-9c93-4992-bfec-f106d8aa0bf6
# ╠═56d88bb9-3b33-4b1a-88ae-d90af4de2bd1
# ╠═c02d0705-6647-4a44-8ae8-fc256f18c4ce
# ╠═d15726ab-5a28-4a24-b5ed-b3c8ecb6c581
# ╠═ae686da9-45d5-4fc2-9cbd-2d828d792407
# ╠═71d3f92e-d143-40dc-8701-37f9053766ef
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
