//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
// Copyright (C) 2010-2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-


#ifndef GRINS_PRIMITIVE_FLOW_FE_VARIABLES_H
#define GRINS_PRIMITIVE_FLOW_FE_VARIABLES_H

// GRINS
#include "grins/fe_variables_base.h"
#include "grins/primitive_flow_variables.h"

namespace GRINS
{

  class PrimitiveFlowFEVariables : public FEVariablesBase,
                                   public PrimitiveFlowVariables
  {
  public:

    PrimitiveFlowFEVariables( const GetPot& input, const std::string& physics_name );
    ~PrimitiveFlowFEVariables(){};

    virtual void init( libMesh::FEMSystem* system );

  protected:

    unsigned int _u_fe_idx, _p_fe_idx;

  private:

    PrimitiveFlowFEVariables();

  };

} // end namespace GRINS

#endif //GRINS_PRIMITIVE_FLOW_VARIABLES_H