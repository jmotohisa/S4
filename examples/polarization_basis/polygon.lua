S = S4.NewSimulation()
--S:SetLattice({1,0}, {-0.5,0.5*math.sqrt(3)})
S:SetLattice({1,0}, {0,1})
--S:SetNumG(160)
S:SetNumG(17)
S:AddMaterial("Silicon", {12,0}) -- real and imag parts
S:AddMaterial("Vacuum", {1,0})

S:AddLayer('AirAbove', 0 , 'Vacuum')
S:AddLayer('Slab', 0.5, 'Silicon')
S:SetLayerPatternPolygon('Slab',     -- which layer to alter
                         'Vacuum',   -- material in polygon
	                     {0.2,0},    -- center
	                     0,          -- tilt angle (degrees)
	                     {0.3, 0.2,  -- vertices
	                      0.3, 0.3,
	                     -0.3, 0.3,
	                     -0.3,-0.3,
	                      0.3,-0.3,
	                      0.3,-0.2,
	                     -0.2,-0.2,
	                     -0.2, 0.2})
S:AddLayerCopy('AirBelow', 0, 'AirAbove')

S:SetExcitationPlanewave(
	{0,0}, -- incidence angles
	{1,0}, -- s-polarization amplitude and phase (in degrees)
	{0,0}) -- p-polarization amplitude and phase

S:SetResolution(16)
S:UsePolarizationDecomposition()
S:SetBasisFieldDumpPrefix('sqr-polygon')

S:SetFrequency(0.4)
S:GetPoyntingFlux('AirBelow', 0)

-- S:OutputLayerPatternRealization('Slab', 64, 64)
