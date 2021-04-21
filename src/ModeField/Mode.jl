# Defines a mode according the horizontal k-vectors and direction in k_z

# mutable struct Mode
#     # Horizontal components of k-vector
#     kXY::_2VectorFloat
# 
#     # Vaccum wavenumber
#     wavenumber::Wavenumber
# 
#     # True is the z-component of the k-vector in the positive z direction.  False otherwise
#     kzPositive::Bool
# 
#     # Mode amplitudes (corresponding to s- and p- modes for propagating waves in unpatterned layers).  Like a Jones vector.  [s,p]
#     A::_2VectorComplex
# 
#     function Mode(kXY::_2VectorFloat, wavenumber::Wavenumber, kzPositive::Bool, A::_2VectorComplex)
#         return new(kXY, wavenumber, kzPositive, A)
#     end
# 
# end
# Mode(k::TU2VectorReal, wavenumber::Wavenumber, kzPositive::Bool, A::TU2VectorComplex) = Mode(_2VectorFloat(k), wavenumber, kzPositive, _2VectorComplex(A))

# Define mode according to real 3-vector for k and a Jones vector polarization
# function Mode(k::_3VectorFloat, A::_2VectorComplex)
# function Mode(k::_3VectorFloat, wavenumber::Wavenumber, A::_2VectorComplex)
#     @assert k[Z] != 0
#     return Mode(k[X:Y], wavenumber, k[Z]>0, A)
# end
# Mode(k::TU3VectorReal, wavenumber::Wavenumber, A::TU2VectorComplex) = Mode(_3VectorFloat(k), wavenumber, _2VectorComplex(A))


# Returns z-component of k-vector of mode
# n is the RI of the layer
# function getkz(mode::Mode, n::Number)::ComplexF64
#     nk₀ = getk₀(mode.wavenumber)*n
#     return sqrt(nk₀^2 - mode.kXY[X]^2 - mode.kXY[Y]^2 ) *bool2posNeg(mode.kzPositive)
# end


# Returns the [kx,ky,kz] vector for this mode
# function kXYtokXYZ(mode::Mode, n::Number)::_3VectorComplex
#     return kXYtokXYZ(mode.kXY, n::Number, mode.wavenumber, mode.kzPositive )
# end

# Get polarization vector
# 𝐏 = Ps*ŝ + Pp*p̂
# function fieldSPtoFieldXYZ(mode::Mode, n::Number)::_3VectorComplex
#     return fieldSPtoFieldXYZ(kXYtokXYZ(mode, n), mode.A)
# end

# function Base.isapprox(mode1::Mode, mode2::Mode; atol::Real=0, rtol::Real=sqrt(eps(Float64)),
#                   nans::Bool=false, norm::Function=abs)::Bool
#     return isapprox(mode1.kXY[X],mode2.kXY[X]; rtol=rtol, atol=atol, nans=nans) &&
#         isapprox(mode1.kXY[Y],mode2.kXY[Y]; rtol=rtol, atol=atol, nans=nans) &&
#         isapprox(mode1.wavenumber.k₀,mode2.wavenumber.k₀; rtol=rtol, atol=atol, nans=nans) &&
#         (mode1.kzPositive == mode1.kzPositive) &&
#         isapprox(mode1.A[X],mode2.A[X]; rtol=rtol, atol=atol, nans=nans) &&
#         isapprox(mode1.A[Y],mode2.A[Y]; rtol=rtol, atol=atol, nans=nans)
# end
