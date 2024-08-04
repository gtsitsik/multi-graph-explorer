classdef graph_tree_node < graph_tree_root
% TODO: Make the slice related properties only available for the first children of the tree.
    properties
        Parent = graph_tree_root.empty
        next_sibling = graph_tree_node.empty
        prev_sibling = graph_tree_node.empty
        size
        type = ''
        is_symmetric = false
        loc
        offset = 0
        is_bipartite = false
        noise_level = 0.0
        sparsity_level = 0.0
        has_global_noise = true
        has_global_sparsity = true
        slices_num = []
        slices_common_entries_perc = 1
        slices_unique_entries_type = ''
    end
    
    
    
    methods
        function delete(node)
            
            % Remove a node from tree.
            if ~isscalar(node)
                error('Input must be scalar')
            end
            prevNode = node.prev_sibling;
            nextNode = node.next_sibling;
            if ~isempty(prevNode)
                prevNode.next_sibling = nextNode;
            end
            if ~isempty(nextNode)
                nextNode.prev_sibling = prevNode;
            end
            node.next_sibling = graph_tree_node.empty;
            node.prev_sibling = graph_tree_node.empty;
            if isvalid(node.Parent)
                node.Parent.Children(node.Parent.Children==node)=[];
            end
            
        end
        
        function insert(node,type,size)
            if ~isempty(node.Parent)
                new_node = graph_tree_node;
                new_node.type = type;
                new_node.size = size;
                
                if ~isempty(node.next_sibling)
                    node.next_sibling.prev_sibling = new_node;
                end
                node.next_sibling = new_node;
                ind = find(node.Parent.Children==node);
                node.Parent.Children = [node.Parent.Children(1:ind) new_node node.Parent.Children(ind+1:end)];
            end
        end
    end
    
    %--------- default examples ----------
    %
    %    methods
    %       function node = dlnode(Data)
    %          % Construct a dlnode object
    %          if nargin > 0
    %             node.Data = Data;
    %          end
    %       end
    %
    %       function insertAfter(newNode, nodeBefore)
    %          % Insert newNode after nodeBefore.
    %          removeNode(newNode);
    %          newNode.Next = nodeBefore.Next;
    %          newNode.Prev = nodeBefore;
    %          if ~isempty(nodeBefore.Next)
    %             nodeBefore.Next.Prev = newNode;
    %          end
    %          nodeBefore.Next = newNode;
    %       end
    %
    %       function insertBefore(newNode, nodeAfter)
    %          % Insert newNode before nodeAfter.
    %          removeNode(newNode);
    %          newNode.Next = nodeAfter;
    %          newNode.Prev = nodeAfter.Prev;
    %          if ~isempty(nodeAfter.Prev)
    %             nodeAfter.Prev.Next = newNode;
    %          end
    %          nodeAfter.Prev = newNode;
    %       end
    %
    %       function removeNode(node)
    %          % Remove a node from a linked list.
    %          if ~isscalar(node)
    %             error('Input must be scalar')
    %          end
    %          prevNode = node.Prev;
    %          nextNode = node.Next;
    %          if ~isempty(prevNode)
    %             prevNode.Next = nextNode;
    %          end
    %          if ~isempty(nextNode)
    %             nextNode.Prev = prevNode;
    %          end
    %          node.Next = dlnode.empty;
    %          node.Prev = dlnode.empty;
    %       end
    %
    %       function clearList(node)
    %          % Clear the list before
    %          % clearing list variable
    %          prev = node.Prev;
    %          next = node.Next;
    %          removeNode(node)
    %          while ~isempty(next)
    %             node = next;
    %             next = node.Next;
    %             removeNode(node);
    %          end
    %          while ~isempty(prev)
    %             node = prev;
    %             prev = node.Prev;
    %             removeNode(node)
    %          end
    %       end
    %    end
    %
    %    methods (Access = private)
    %       function delete(node)
    %          clearList(node)
    %       end
    %    end
end