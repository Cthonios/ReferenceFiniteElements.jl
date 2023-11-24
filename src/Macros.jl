abstract type RefFEMacroException <: Exception end

struct ValuesException{T <: Type, A <: AbstractArray} <: RefFEMacroException
  type::T
  array::A
end

Base.show(io::IO, e::ValuesException) = println(io, 
  "Error initializing array of type $(e.type) for array $(e.array)."
)

values_error(type, array) = throw(ValuesException(type, array))

macro ValuesArray(type, array)
  if type == :SVector
    :(
      @SVector $array
    )
  elseif type == :MVector
    :(
      @MVector $array
    )
  elseif type == :Vector
    :(
      $array
    )
  else
    values_error(type, array)
  end
end