module AdaptExt

import Adapt
import Adapt: adapt
import ReferenceFiniteElements as RFE

function Adapt.adapt_structure(to, re::RFE.ReferenceFE)
    bdofs = adapt(to, re.boundary_dofs)
    bnorms = adapt(to, re.boundary_normals)
    cell_interps = adapt(to, re.cell_interps)
    surf_interps = adapt(to, re.surf_interps)
    return RFE.ReferenceFE{typeof(re.element), RFE.quadrature_rule(re)}(
        re.element, bdofs, bnorms, cell_interps, surf_interps
    )
end

end # module
