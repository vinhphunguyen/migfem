classdef patch2D < handle
    % This class implements a 2D NURBS patch.
    % Used in a multi-patches IGA code.
    %
    % Vinh Phu Nguyen, February 2012
    % nvinhphu@gmail.com
    
    properties
        % geometry related
        
        uKnot
        vKnot
        p
        q
        controlPts
        weights
        noPtsX
        noPtsY
        
        % analysis related
        
        element
        elementLocal 
        index
        elRangeU
        elRangeV
        elConnU
        elConnV
        noElemsU
        noElemsV
    end
    
    methods
        % class constructor
        function obj = patch2D(uKnot,vKnot,p,q,pts,weights)
            if(nargin > 0)
                obj.uKnot      = uKnot;
                obj.vKnot      = vKnot;
                obj.p          = p;
                obj.q          = q;
                obj.controlPts = pts;
                obj.weights    = weights;
                
                uniqueUKnots   = unique(uKnot);
                uniqueVKnots   = unique(vKnot);                
                obj.noElemsU   = length(uniqueUKnots)-1; 
                obj.noElemsV   = length(uniqueVKnots)-1;  
                obj.noPtsX     = length(uKnot)-p-1;
                obj.noPtsY     = length(vKnot)-q-1;
            end
        end
        
        % mesh generation
        ret = generateMesh (obj, node_pattern, node_patternLocal)
    end
end