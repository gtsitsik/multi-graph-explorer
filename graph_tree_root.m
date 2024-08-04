classdef graph_tree_root < handle
    properties
        Children = graph_tree_node.empty
        Data = [];
        labels = [];
    end
    
    methods
        
        function G_out = deep_copy(G_in)
            props = properties(G_in);
            G_out = feval(metaclass(G_in).Name);
            for cur_prop = props'
                switch cur_prop{1}
                    case 'Children'
                        for i = 1:numel(G_in.Children)
                            G_out.Children(i) = G_in.Children(i).deep_copy();
                        end
                    case {'Parent','next_sibling','prev_sibling'}
                        continue;
                    otherwise
                        G_out.(cur_prop{1}) = G_in.(cur_prop{1});
                end
            end
        end
        
        
        function bool = is_equivalent_to(G_in,G_out)
            props = properties(G_in);
            bool = true;
            for cur_prop = props'
                switch cur_prop{1}
                    case 'Children'
                        if numel(G_in.Children)~=numel(G_out.Children)
                            bool = false;
                            return
                        end
                        for i = 1:numel(G_in.Children)
                            if ~is_equivalent_to(G_in.Children(i),G_out.Children(i))
                                bool = false;
                                return
                            end
                        end
                    case {'Parent','next_sibling','prev_sibling'}
                        continue;
                    otherwise
                        if ~isequal(G_out.(cur_prop{1}),G_in.(cur_prop{1}))
                            bool = false;
                            return
                        end
                end
            end
        end
        
    end
end