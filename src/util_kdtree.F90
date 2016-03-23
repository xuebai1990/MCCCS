MODULE util_kdtree
    use var_type, only:dp
    use util_runtime,only:err_exit
    use util_random,only:random
    use sim_system,only:interval,tree_node,tree
    implicit none
    private

    ! public subroutines
    public :: init_tree, insert_node, range_search, find_min, find_max, delete_node, &
             check_tree_coord, check_tree_cube, search_node_in_tree, update_cube, &
             init_tree_coordinates, update_tree_height, vector_dist_square, &
             empty_tree, allocate_kdtree, construct_kdtree, update_box_kdtree, read_kdtree

    ! global variables
    real, allocatable, public :: search_output(:,:) !< the returned array in the range search 
    integer, public :: current_dim !< the current dimension of the search_output
    real, public :: rmin_square, rcut_square
    real :: coord_replaced(3) !< used in the deletion 

contains

    ! Calculate the distance squared between two vectors
    ! vector1, vector2: two vectors being calculated
    ! max_dist_square: if the distance exceeds this value, no need to further calculate
    ! dist_square: return the distance if it does not exceed the max_dist_square
    function vector_dist_square(vector1, vector2, max_dist_square) result(dist_square)
        real, dimension(3), intent(in) :: vector1, vector2
        real :: max_dist_square, dist_square, vector
        integer :: i

        dist_square = 0.0E0_dp
        do i = 1, 3
            vector = vector1(i) - vector2(i) 
            dist_square = dist_square + vector * vector
            if (dist_square .gt. max_dist_square) return
        end do
    end function vector_dist_square

    ! Calculate the distance squared between one point and the cube of a node
    ! mol_tree: the pointer the tree
    ! vector: the coordinate
    ! current_node: the pointer to the node whose cube is being used in the calculation
    ! max_dist_square: if the distance exceeds this value, no need to further calculate
    ! dist_square: return the distance if it does not exceed the max_dist_square
    function dist_square_from_point_to_cube(mol_tree, vector, current_node, max_dist_square) result(dist_square)
        type(tree), pointer :: mol_tree
        type(tree_node), pointer :: current_node 
        real, dimension(3), intent(in) :: vector
        real :: dist_square, max_dist_square, upper, lower
        integer :: i_dim, cut_dim

        dist_square = 0.0E0_dp
        do i_dim = 1, 3
            if (i_dim .eq. 1) then
                cut_dim = cut_dim_prev(current_node%cut_dim)
                lower = current_node%cube%lower
                upper = current_node%cube%upper
            else if (i_dim .eq. 2) then
                cut_dim = cut_dim_prev(cut_dim) 
                if (associated(current_node%parent_node)) then
                    lower = current_node%parent_node%cube%lower
                    upper = current_node%parent_node%cube%upper
                else
                    lower = mol_tree%cube(cut_dim)%lower
                    upper = mol_tree%cube(cut_dim)%upper
                end if
            else
                cut_dim = cut_dim_prev(cut_dim)
                if (current_node%height .ge. 3) then
                    lower = current_node%parent_node%parent_node%cube%lower
                    upper = current_node%parent_node%parent_node%cube%upper
                else
                    lower = mol_tree%cube(cut_dim)%lower
                    upper = mol_tree%cube(cut_dim)%upper
                end if
            end if
    
            if (vector(cut_dim) .lt. lower) then
                dist_square = dist_square + (vector(cut_dim) - lower)**2
            else if (vector(cut_dim) .gt. upper) then
                dist_square = dist_square + (vector(cut_dim) - upper)**2
            end if
            
            if (dist_square .gt. max_dist_square) return
        end do

    end function

    ! Calculate the next dimension
    ! cut_dim: the cutting dimension
    ! next_cut_dim: the next cutting dimension
    function cut_dim_next(cut_dim) result(next_cut_dim)
        integer, intent(in) :: cut_dim
        integer :: next_cut_dim
        
        if (cut_dim .eq. 1) then
            next_cut_dim = 2
        else if (cut_dim .eq. 2) then
            next_cut_dim = 3
        else
            next_cut_dim = 1
        end if
    end function

    ! Return the previous dimension
    ! cut_dim: the cutting dimension
    ! previous_cut_dim: the previous cutting dimension
    function cut_dim_prev(cut_dim) result(previous_cut_dim)
        integer, intent(in) :: cut_dim
        integer :: previous_cut_dim

        if (cut_dim .eq. 1) then
            previous_cut_dim = 3
        else if (cut_dim .eq. 2) then
            previous_cut_dim = 1
        else
            previous_cut_dim = 2
        end if
    end function cut_dim_prev
        
    ! Determine whether to go left or right from a parent node
    ! coord_to_add: 3d real array, the coordinate of the node to be added/searched/deleted
    ! coord_to_compare: 3d real array, the coodinate of the node to be compared (typically the current node)
    ! cut_dim: the cutting dimension of the current node
    ! lLeft: return true if it should go left, false if it should go right
    function l_left(coord_to_add, coord_to_compare, cut_dim) result(lLeft)
        real, dimension(3) :: coord_to_add, coord_to_compare
        integer :: cut_dim, next_cut_dim, next_next_cut_dim
        logical :: lLeft

        if (coord_to_add(cut_dim) .lt. coord_to_compare(cut_dim)) then
            lLeft = .true.
        else if (coord_to_add(cut_dim) .gt. coord_to_compare(cut_dim)) then
            lLeft = .false.
        else
            next_cut_dim = cut_dim_next(cut_dim)
            if (coord_to_add(next_cut_dim) .lt. coord_to_compare(next_cut_dim)) then
                lLeft = .true.
            else if (coord_to_add(next_cut_dim) .gt. coord_to_compare(next_cut_dim)) then
                lLeft = .false.
            else
                next_next_cut_dim = cut_dim_next(next_cut_dim)
                if (coord_to_add(next_next_cut_dim) .lt. coord_to_compare(next_next_cut_dim)) then
                    lLeft = .true.
                else
                    lLeft = .false.
                end if
            end if
        end if
    end function l_left


    ! Set the cube for current_node
    ! mol_tree: the pointer to the tree
    ! current_node: the pointer to the current node
    ! lLeft: whether the current node is the left or right child of its parent
    subroutine set_cube(mol_tree, current_node, lLeft)
        type(tree), pointer :: mol_tree
        type(tree_node), pointer :: current_node
        integer :: cut_dim
        logical, intent(in) :: lLeft
        real :: lower, upper

        current_node%l_cube_updated = .true.

        ! if it's the root, the cube is for z-dimension
        if (.not. associated(current_node%parent_node)) then
            current_node%cube%lower = mol_tree%cube(3)%lower
            current_node%cube%upper = mol_tree%cube(3)%upper
            return
        end if

        cut_dim = current_node%parent_node%cut_dim !< cut_dim is the parent's cut_dim

        ! Determine whether there is great grandparent for the current node
        if (current_node%height .ge. 4) then
            lower = current_node%parent_node%parent_node%parent_node%cube%lower
            upper = current_node%parent_node%parent_node%parent_node%cube%upper
        else
            lower = mol_tree%cube(cut_dim)%lower
            upper = mol_tree%cube(cut_dim)%upper
        end if

        if (lLeft) then
            current_node%cube%lower = lower
            current_node%cube%upper = current_node%parent_node%coord(cut_dim)
        else    
            current_node%cube%lower = current_node%parent_node%coord(cut_dim)
            current_node%cube%upper = upper 
        end if
        
    end subroutine set_cube

    ! to initialize a kd-tree
    ! mol_tree: the pointer to the tree being calculated
    ! ibox: the box that this tree represents
    ! cubeP: the cube of the root
    ! mol_tree_initialized: return the pointer to the tree being initialized
    function init_tree(mol_tree, ibox, cubeP) result(mol_tree_initialized)
        type(tree), pointer :: mol_tree, mol_tree_initialized
        integer :: ibox
        type(interval), pointer :: cubeP(:)
        if (associated(mol_tree)) call err_exit(__FILE__,__LINE__,'kd-tree already existed!',1)
        allocate(mol_tree)
        allocate(mol_tree%cube(3))
        allocate(mol_tree%bound(3))
        mol_tree%cube = cubeP
        mol_tree%node_num = 0
        mol_tree%height = 0
        mol_tree%box = ibox
        nullify(mol_tree%tree_root)
        mol_tree_initialized => mol_tree
    end function init_tree

    ! empty all the nodes in the tree
    subroutine empty_tree(mol_tree)
        type(tree), pointer :: mol_tree
        if (associated(mol_tree)) call empty_node(mol_tree%tree_root)
        deallocate(mol_tree%cube)
        deallocate(mol_tree%bound)
        deallocate(mol_tree)
    end subroutine empty_tree

    ! recursively empty the nodes
    recursive subroutine empty_node(current_node)
        type(tree_node), pointer :: current_node
        if (associated(current_node%left_node)) then
            call empty_node(current_node%left_node)
            nullify(current_node%left_node)
        end if

        if (associated(current_node%right_node)) then
            call empty_node(current_node%right_node)
            nullify(current_node%right_node)
        end if

        if (associated(current_node%cube)) deallocate(current_node%cube) 
        deallocate(current_node) 
    end subroutine empty_node
   
    ! Insert a node to the tree
    ! mol_tree: the pointer to the tree where the node is being inserted
    ! coord_to_add: 3d real array, the coordinate of the node to be added
    ! ichain: the chain number of the molecule
    ! ibead: the bead number of the molecule 
    ! ix, iy, iz: which periodic image it is on the x, y or z axis, a value between -1, 1 and 0
    ! mol_tree_inserted: return the pointer to the tree after the insertion
    function insert_node(mol_tree, coord_to_add, ichain, ibead, ix, iy, iz) result(mol_tree_inserted)
        type(tree), pointer :: mol_tree, mol_tree_inserted
        real, dimension(3) :: coord_to_add
        integer, intent(in) :: ichain, ibead, ix, iy, iz
    
        if (.not. associated(mol_tree)) call err_exit(__FILE__,__LINE__,'kd-tree has not been initialized yet!',1)

        mol_tree%tree_root => insert_node_into_tree(mol_tree, coord_to_add, ichain, ibead, ix, iy, iz, &
                                mol_tree%tree_root, 1, 1, .true., null())
        mol_tree_inserted => mol_tree
    end function insert_node

    ! Insert the coordinates coord_to_add to the proper node
    ! mol_tree: pointer to the tree
    ! coord_to_add: 3d real array, the coordinate of the node to be added
    ! ichain: the chain number of the molecule
    ! ibead: the bead number of the molecule
    ! ix, iy, iz: which periodic image it is on the x, y or z axis, a value between -1, 1 and 0 
    ! current_node: the current searching node
    ! height: the height of the current node
    ! cut_dim: the cutting dimension of the current node
    ! lLeft: true if it's the left of its parent node
    ! parent_node: the pointer to the parent node
    ! node_added: return the current searching node
    recursive function insert_node_into_tree(mol_tree, coord_to_add, ichain, ibead, ix, iy, iz, current_node, height&
                        ,cut_dim, lLeft, parent_node) result(node_added)
        type(tree), pointer :: mol_tree
        real, dimension(3) :: coord_to_add
        integer :: ichain, ibead, ix, iy, iz
        type(tree_node), pointer :: current_node, parent_node, node_added
        integer :: height
        integer, intent(in) :: cut_dim
        logical, intent(in) :: lLeft
        integer :: next_cut_dim
        next_cut_dim = cut_dim_next(cut_dim)

        if (.not. associated(current_node)) then
            allocate(current_node)
            current_node%coord = coord_to_add
            current_node%ichain = ichain
            current_node%ibead = ibead
            current_node%ix = ix
            current_node%iy = iy
            current_node%iz = iz
            current_node%l_cube_updated = .false.
            current_node%parent_node => parent_node
            current_node%cut_dim = cut_dim
            current_node%height = height
            !current_node%cube%lower = 0
            !current_node%cube%upper = 0
            allocate(current_node%cube)
            if (.not. associated(current_node%parent_node)) call set_cube(mol_tree, current_node, lLeft)
            mol_tree%node_num = mol_tree%node_num + 1
            if (height .gt. mol_tree%height) mol_tree%height = height
            nullify(current_node%left_node, current_node%right_node)
        else if ((vector_dist_square(current_node%coord, coord_to_add, 1e-6) .lt. 1e-6) &
                .and. (current_node%ichain .eq. ichain) .and. (current_node%ibead .eq. ibead) &
                .and. (current_node%ix .eq. ix) .and. (current_node%iy .eq. iy) .and. (current_node%iz .eq. iz)) then
            write(*,*) "current_node: coordinates", current_node%coord
            write(*,*) "current_node: ichain, ibead", current_node%ichain, current_node%ibead
            write(*,*) "current_node: ix, iy, iz", current_node%ix, current_node%iy, current_node%iz
            write(*,*) "coord_to_add: coordinates", coord_to_add
            write(*,*) "coord_to_add: ichain, ibead", ichain, ibead
            write(*,*) "coord_to_add: ix, iy, iz", ix, iy, iz
            call err_exit(__FILE__,__LINE__,'Duplicate keys added when constructing the kd-tree',1)
        else if (l_left(coord_to_add, current_node%coord, cut_dim)) then 
            current_node%left_node => insert_node_into_tree(mol_tree, coord_to_add, ichain, ibead &
                , ix, iy, iz, current_node%left_node, height+1, next_cut_dim, .true., current_node)
        else
            current_node%right_node => insert_node_into_tree(mol_tree, coord_to_add, ichain, ibead &
                , ix, iy, iz, current_node%right_node, height+1, next_cut_dim, .false., current_node)
        end if

        node_added => current_node

    end function insert_node_into_tree
    
    ! Return the minimum node in the tree at dimension min_dim
    ! current_node: the current node being examined
    ! min_dim: the dimension where the minimum is
    ! cut_dim: the current dimension
    ! min_node: return the pointer to the minimum node
    recursive function find_min(current_node, min_dim, cut_dim) result(min_node)
        type(tree_node), pointer :: current_node, min_node, min_nodeLeft, min_node_right
        integer, intent(in) :: min_dim, cut_dim
        real :: min_node_coord(3)
        logical :: lLeft, lRight
        
        if (.not. associated(current_node)) min_node => null()
        
        if (cut_dim .eq. min_dim) then
            if (.not. associated(current_node%left_node)) then 
                min_node => current_node
            else 
                min_node => find_min(current_node%left_node, min_dim, cut_dim_next(cut_dim))
            end if
        else
            min_node => current_node
            min_node_coord = current_node%coord
            lLeft = associated(current_node%left_node)
            lRight = associated(current_node%right_node)

            ! if there is sub-tree, compare the value with left and right sub-tree minimums
            if (lLeft) then
                min_nodeLeft => find_min(current_node%left_node, min_dim, cut_dim_next(cut_dim))
                if (l_left(min_nodeLeft%coord, min_node_coord, min_dim)) then
                    min_node => min_nodeLeft
                    min_node_coord = min_node%coord
                end if
            end if

            if (lRight) then
                min_node_right => find_min(current_node%right_node, min_dim, cut_dim_next(cut_dim))
                if (l_left(min_node_right%coord, min_node_coord, min_dim)) min_node => min_node_right
            end if
        end if

    end function find_min

    ! Return the maximum node in the tree at dimension min_dim
    ! current_node: the current node being examined
    ! max_dim: the dimension where the maximum is
    ! cut_dim: the current dimension
    ! max_node: return the pointer to the maximum node
    recursive function find_max(current_node, max_dim, cut_dim) result(max_node)
        type(tree_node), pointer :: current_node, max_node, max_node_left, max_nodeRight
        integer, intent(in) :: max_dim, cut_dim
        real :: max_nodeCoord(3)
        logical :: lLeft, lRight
        
        if (.not. associated(current_node)) max_node => null()
        
        if (cut_dim .eq. max_dim) then
            if (.not. associated(current_node%right_node)) then 
                max_node => current_node
            else 
                max_node => find_max(current_node%right_node, max_dim, cut_dim_next(cut_dim))
            end if
        else
            max_node => current_node
            max_nodeCoord = current_node%coord
            lLeft = associated(current_node%left_node)
            lRight = associated(current_node%right_node)

            ! if there is sub-tree, compare the value with left and right sub-tree minimums
            if (lLeft) then
                max_node_left => find_max(current_node%left_node, max_dim, cut_dim_next(cut_dim))
                if (.not. l_left(max_node_left%coord, max_nodeCoord, max_dim)) then
                    max_node => max_node_left
                    max_nodeCoord = max_node%coord
                end if
            end if

            if (lRight) then
                max_nodeRight => find_max(current_node%right_node, max_dim, cut_dim_next(cut_dim))
                if (.not. l_left(max_nodeRight%coord, max_nodeCoord, max_dim)) max_node => max_nodeRight
            end if
        end if

    end function find_max

    ! Update the height of the tree
    ! mol_tree: the pointer to the tree
    subroutine update_tree_height(mol_tree)
        type(tree), pointer :: mol_tree
    
        mol_tree%height = 0
        if (.not. associated(mol_tree)) return
        call update_height_in_tree(mol_tree, mol_tree%tree_root, 1, 1)
    end subroutine update_tree_height

    ! Recursively find the height of the tree
    ! mol_tree: the pointer to the tree
    ! current_node: the current node being worked on
    ! cut_dim: the cutting dimension of the current node
    ! height: the height of the current node
    recursive subroutine update_height_in_tree(mol_tree, current_node, cut_dim, height)
        type(tree), pointer :: mol_tree
        type(tree_node), pointer :: current_node
        integer :: cut_dim, next_cut_dim, height

        next_cut_dim = cut_dim_next(cut_dim)

        if (height .gt. mol_tree%height) mol_tree%height = height
        if (associated(current_node%left_node)) call update_height_in_tree(mol_tree, current_node%left_node, next_cut_dim, height+1)
        if (associated(current_node%right_node)) call update_height_in_tree(mol_tree, current_node%right_node, next_cut_dim, height+1)
         
    end subroutine update_height_in_tree

    ! Delete a node in the tree, return the pointer to the tree
    ! mol_tree: the pointer to the tree
    ! coord_to_delete: 3d coordinates of the nodes to be deleted
    ! mol_tree_deleted: return the pointer to the tree
    function delete_node(mol_tree, coord_to_delete) result(mol_tree_deleted)
        type(tree), pointer :: mol_tree, mol_tree_deleted
        real, dimension(3), intent(in) :: coord_to_delete
        real, dimension(3) :: coord_to_compare
        type(tree_node), pointer :: node_replaced
        !type(tree_node), pointer :: node_to_delete
        !type(tree_node), pointer, optional :: node_to_delete_input
        !logical :: l_node_to_delete

        !if (present(node_to_delete_input)) then
        !    l_node_to_delete = .true.
        !    node_to_delete => node_to_delete_input
        !else
        !    l_node_to_delete = .false.
        !end if

        coord_replaced = [10000, 10000, 10000] !< set to unrealistic values       
        coord_to_compare = coord_replaced

        !if (l_node_to_delete) then
        !   ! if the node to be deleted has been specified
        !    node_to_delete => delete_node_in_tree(mol_tree, node_to_delete, coord_to_delete, node_to_delete%cut_dim, .false.)
        !else
            ! if not, start searching from the tree root
        mol_tree%tree_root => delete_node_in_tree(mol_tree, mol_tree%tree_root, coord_to_delete, 1, .false.)      
        !end if
 
        mol_tree%node_num = mol_tree%node_num - 1
        mol_tree_deleted => mol_tree
      
        if (vector_dist_square(coord_to_compare, coord_replaced, 1e-6) .gt. 1e-6) then
            node_replaced => search_node_in_tree(mol_tree, mol_tree%tree_root, coord_replaced, 1)
            node_replaced => reset_cube_status(node_replaced)
        end if

    end function delete_node

    ! Recursively delete a node in the tree
    ! mol_tree: the pointer to the tree
    ! current_node: the node currently being examined
    ! coord_to_delete: 3d coordinates of the nodes to be deleted
    ! cut_dim: the cutting dimension of the current node
    ! l_update_cube: whether to update cube, false until the the first occurence of node deletion
    ! node_deleted: return the pointer to the node AFTER the deletion
    recursive function delete_node_in_tree(mol_tree, current_node, coord_to_delete, cut_dim, l_update_cube) result(node_deleted)
        type(tree), pointer :: mol_tree
        type(tree_node), pointer :: current_node, node_deleted, node_replaced
        real, dimension(3), intent(in) :: coord_to_delete
        integer, intent(in) :: cut_dim
        logical, intent(in) :: l_update_cube
        integer :: next_cut_dim

        next_cut_dim = cut_dim_next(cut_dim)

        if (.not. associated(current_node)) then
            write(*,*) coord_to_delete
            call err_exit(__FILE__,__LINE__,'Did not find the node!',1)
        end if

        if (vector_dist_square(current_node%coord, coord_to_delete, 1e-6) .lt. 1e-6) then

            ! if this is the node to delete
            if (associated(current_node%right_node)) then
                ! if there's only right node, use the minimum from the right sub-tree
                node_replaced => find_min(current_node%right_node, cut_dim, next_cut_dim)

                current_node%coord = node_replaced%coord
                current_node%ichain = node_replaced%ichain
                current_node%ibead = node_replaced%ibead
                current_node%ix = node_replaced%ix
                current_node%iy = node_replaced%iy
                current_node%iz = node_replaced%iz
                
                ! record the coordinates to coord_replaced if l_update_cube == .false.
                if (.not. l_update_cube) coord_replaced = node_replaced%coord

                current_node%right_node => delete_node_in_tree(mol_tree, current_node%right_node, current_node%coord, &
                                        next_cut_dim, .true.)
            else if (associated(current_node%left_node)) then
                ! if there's only left node, swap sub-trees, and use the minimum from the next right-tree
                node_replaced => find_min(current_node%left_node, cut_dim, next_cut_dim)
                current_node%coord = node_replaced%coord
                current_node%ichain = node_replaced%ichain
                current_node%ibead = node_replaced%ibead
                current_node%ix = node_replaced%ix
                current_node%iy = node_replaced%iy
                current_node%iz = node_replaced%iz
                
                ! record the coordinates to coord_replaced if l_update_cube == .false.
                if (.not. l_update_cube) coord_replaced = node_replaced%coord

                current_node%right_node => delete_node_in_tree(mol_tree, current_node%left_node, current_node%coord, &
                                        next_cut_dim, .true.)
                current_node%left_node => null()
            else
                ! if it's only a leaf, just remove
                deallocate(current_node%cube)
                deallocate(current_node)
                node_deleted => null()
                return
            end if

        ! if this is not the node to delete, search for it
        else if (l_left(coord_to_delete,  current_node%coord, cut_dim)) then
            current_node%left_node => delete_node_in_tree(mol_tree, current_node%left_node, coord_to_delete, next_cut_dim, l_update_cube)
        else
            current_node%right_node => delete_node_in_tree(mol_tree, current_node%right_node, coord_to_delete, next_cut_dim, l_update_cube)
        end if

        ! Return the node 
        node_deleted => current_node

    end function delete_node_in_tree

    ! Reset all the l_cube_updated to be .false. recursively
    ! current_node: the starting node to be updated
    ! node_updated: return the updated node pointer
    recursive function reset_cube_status(current_node) result(node_updated)
        type(tree_node), pointer :: current_node, node_updated

        current_node%l_cube_updated = .false.
        if (associated(current_node%left_node)) current_node%left_node => reset_cube_status(current_node%left_node)        
        if (associated(current_node%right_node)) current_node%right_node => reset_cube_status(current_node%right_node)        

        node_updated => current_node
    end function reset_cube_status

    ! Update the cube for all the sub-trees of coord_to_update
    ! mol_tree: the pointer to the tree
    ! coord_to_update: 3d coordinates of the node whose sub-trees need to be updated
    ! mol_treeUpdated: return the pointer to the tree after the cube update
    function update_cube(mol_tree, coord_to_update) result(mol_treeUpdated)
        type(tree), pointer :: mol_tree, mol_treeUpdated
        real, dimension(3) :: coord_to_update
        type(tree_node), pointer :: node_to_update
        
        node_to_update => search_node_in_tree(mol_tree, mol_tree%tree_root, coord_to_update, 1)
        node_to_update => update_cube_in_tree(mol_tree, node_to_update)
        mol_treeUpdated => mol_tree
    end function update_cube

    ! Update the cube for all the sub-trees of coord_to_update
    ! mol_tree: the pointer to the tree
    ! current_node: the node currently working on
    ! node_updated: the updated current node
    recursive function update_cube_in_tree(mol_tree, current_node) result(node_updated)
        type(tree), pointer :: mol_tree
        type(tree_node), pointer :: current_node, node_updated

        if (associated(current_node%left_node)) then
            call set_cube(mol_tree, current_node%left_node, .true.)
            current_node%left_node => update_cube_in_tree(mol_tree, current_node%left_node)
        end if

        if (associated(current_node%right_node)) then
            call set_cube(mol_tree, current_node%right_node, .false.)
            current_node%right_node => update_cube_in_tree(mol_tree, current_node%right_node)
        end if

        node_updated => current_node

    end function update_cube_in_tree

    ! Search a node in the tree
    ! current_node: the node currently working on
    ! coord_to_search: 3d coordinates of the node to be searched
    ! cut_dim: the cutting dimension of the current_node
    ! node_found: the node found
    recursive function search_node_in_tree(mol_tree, current_node, coord_to_search, cut_dim) result(node_found)
        type(tree), pointer :: mol_tree
        type(tree_node), pointer :: current_node, node_found
        real, dimension(3) :: coord_to_search
        integer, intent(in) :: cut_dim

        if (.not. associated(current_node)) then
            node_found => null()
            return
        end if

        if (vector_dist_square(current_node%coord, coord_to_search, 1e-6) .le. 1e-6) then
            node_found => current_node
        else if (l_left(coord_to_search, current_node%coord, cut_dim)) then
            node_found => search_node_in_tree(mol_tree, current_node%left_node, coord_to_search, cut_dim_next(cut_dim))
        else    
            node_found => search_node_in_tree(mol_tree, current_node%right_node, coord_to_search, cut_dim_next(cut_dim))
        end if
        
    end function search_node_in_tree

    ! Check the validity of the tree, for the coordinates only
    ! mol_tree: the pointer to the tree
    ! lValid: whether the tree structure is correct
    function check_tree_coord(mol_tree) result(lValid)
        type(tree), pointer :: mol_tree
        logical :: lValid
        lValid = check_coord_in_tree(mol_tree%tree_root, 1)
    end function check_tree_coord
    
    ! Check the validity of the node in the tree
    ! current_node: the node being examined
    ! lValid: whether the tree structure at this node is correct
    recursive function check_coord_in_tree(current_node, cut_dim) result(lValid)
        type(tree_node), pointer :: current_node, max_node_left, min_node_right
        logical :: lValid
        integer :: cut_dim, next_cut_dim

        ! if the cutting dimension does not match
        if (current_node%cut_dim .ne. cut_dim) then
            lValid = .false.
            return
        end if

        next_cut_dim = cut_dim_next(cut_dim)
        lValid = .true.

        if (associated(current_node%left_node)) then
            max_node_left => find_max(current_node%left_node, cut_dim, next_cut_dim)
            if (l_left(current_node%coord, max_node_left%coord, cut_dim)) then
                lValid = .false.
                return
            else
                lValid = check_coord_in_tree(current_node%left_node, next_cut_dim)
                if (.not. lValid) return
            end if
        end if
    
        if (associated(current_node%right_node)) then
            min_node_right => find_min(current_node%right_node, cut_dim, next_cut_dim)
            if (.not. l_left(current_node%coord, min_node_right%coord, cut_dim)) then
                lValid = .false.
                return
            else
                lValid = check_coord_in_tree(current_node%right_node, next_cut_dim)
                if (.not. lValid) return
            end if
        end if
    end function check_coord_in_tree

    ! Check the validity of the cube in the tree
    ! mol_tree: the pointer to the tree
    ! lValid: whetehr the tree cube is correct
    function check_tree_cube(mol_tree) result(lValid)
        type(tree), pointer :: mol_tree
        logical :: lValid
        lValid = check_cube_in_tree(mol_tree, mol_tree%tree_root)
    end function check_tree_cube

    ! Check the validity of the cube in the current_node
    ! mol_tree: the pointer to the tree
    ! current_node: the node being examined
    ! lValid: whether the tree cube at this node is correct
    recursive function check_cube_in_tree(mol_tree, current_node) result(lValid)
        type(tree), pointer :: mol_tree
        type(tree_node), pointer :: current_node
        logical :: lValid
        integer :: previous_cut_dim
        real :: upper, lower

        ! if it's the root
        if (.not. associated(current_node%parent_node)) then
            if ((current_node%cube%upper .eq. mol_tree%cube(3)%upper) &
                    .and. (current_node%cube%lower .eq. mol_tree%cube(3)%lower)) then
                lValid = .true.
            else
                lValid = .false.
                return
            end if
 
            !if (.not. lValid) return
            if (associated(current_node%left_node)) lValid = check_cube_in_tree(mol_tree, current_node%left_node)
            if (.not. lValid) return
            if (associated(current_node%right_node)) lValid = check_cube_in_tree(mol_tree, current_node%right_node)
            return
        end if

        previous_cut_dim = current_node%parent_node%cut_dim

        ! if it's NOT the root, check the cube
        lower = current_node%cube%lower
        upper = current_node%cube%upper

        if (l_left(current_node%coord, current_node%parent_node%coord, previous_cut_dim)) then
            call set_cube(mol_tree, current_node, .true.)
        else
            call set_cube(mol_tree, current_node, .false.)
        end if

        if ((lower .ne. current_node%cube%lower) .or. (upper .ne. current_node%cube%upper)) then
            lValid = .false.
            return
        else
            lValid = .true.
        end if       

        ! recursively check its leaves
        if (associated(current_node%left_node)) lValid = check_cube_in_tree(mol_tree, current_node%left_node)
        if (.not. lValid) return
        if (associated(current_node%right_node)) lValid = check_cube_in_tree(mol_tree, current_node%right_node)
        return
        
    end function check_cube_in_tree

    ! Search for the coordinates within rcut with respect to a reference point ref_coord
    ! Return loverlap = true if there is an r larger than rmin
    ! otherwise return an array of eligible search_output
    ! mol_tree: pointer the tree
    ! ref_coord: the reference coordinate, whose interactions with all other beads are being calculated
    ! max_dim: maximum dimension of the search_output
    ! rmin: minimum distance allowed
    ! rcut: cutoff distance
    ! ichain: ref_coord is for ichain; do not search for bead in ichain
    ! ibead: ref_coord is for ibead in ichain; used to determine what type of bead it is and whether it has LJ or QQ interaction
    ! lsumup: true if called from sumup suborutine, accout for only interactions when chain number is greater than ichain
    !         false if called from energy or boltz subroutine, account for all the interactions when bead is not in ichain
    ! loverlap: return true if there's an overlap
    ! actual_dim: actual dimension of the search_output
    ! search_output_array: 3*actual_dim real array which have information about the qualified beads to calculate interactions
    ! dist_calc_num: number of distance calculation used
    ! lPressure: whether called from pressure calculation, if so, have additional dimension of output arrays for r*uij
    subroutine range_search(mol_tree, ref_coord, max_dim, rmin, rcut, ichain, ibead, lsumup, &
        loverlap, actual_dim, search_output_array, dist_calc_num, lPressure)
        type(tree), pointer :: mol_tree
        real, dimension(3), intent(in) :: ref_coord
        integer, intent(in) :: max_dim, ichain, ibead 
        real, intent(in) :: rmin, rcut
        logical, intent(in) :: lPressure
        logical :: lsumup, loverlap
        integer :: actual_dim
        real, allocatable :: search_output_array(:,:)
        integer :: i, dist_calc_num, search_output_dim

        ! Initialize the variables
        if (allocated(search_output)) deallocate(search_output)

        !< dist_square, ichain, ibead for 3 dimensions, if lPressure, additional 3d for rxuij, ryuij and rzuij
        if (lPressure) then
            search_output_dim = 6
        else
            search_output_dim = 3
        end if

        allocate(search_output(search_output_dim, max_dim))

        current_dim = 1
        loverlap = .false.
        rmin_square = rmin * rmin
        rcut_square = rcut * rcut
        dist_calc_num = 0
        
        ! Start searching
        call range_search_in_tree(mol_tree, mol_tree%tree_root, ref_coord, 1, ichain, ibead, lsumup, loverlap, dist_calc_num, lPressure)
        actual_dim = current_dim - 1
        if (allocated(search_output_array)) deallocate(search_output_array)

        allocate(search_output_array(search_output_dim, actual_dim))

        do i = 1, actual_dim
            search_output_array(1:search_output_dim, i) = search_output(1:search_output_dim, i)
        end do
    end subroutine range_search

    ! recursively do the range search
    ! mol_tree: pointer the tree
    ! current_node: the node currently working on
    ! ref_coord: the reference coordinate, whose interactions with all other beads are being calculated
    ! cut_dim: the cutting dimension of the current node
    ! ichain: ref_coord is for ichain; do not search for bead in ichain
    ! ibead: ref_coord is for ibead in ichain; used to determine what type of bead it is and whether it has LJ or QQ interaction
    ! lsumup: true if called from sumup suborutine, accout for only interactions when chain number is greater than ichain
    !         false if called from energy or boltz subroutine, account for all the interactions when bead is not in ichain
    ! loverlap: return true if there's an overlap
    ! dist_calc_num: number of distance calculation used
    ! lPressure: whether called from pressure calculations, add additional 3d output for r*uij
    recursive subroutine range_search_in_tree(mol_tree, current_node, ref_coord, cut_dim, ichain, ibead, lsumup&
                , loverlap, dist_calc_num, lPressure)
        use sim_system,only:io_output

        type(tree), pointer :: mol_tree
        real, dimension(3), intent(in) :: ref_coord
        integer, intent(in) :: cut_dim, ichain, ibead
        type(tree_node), pointer :: current_node, node_closer, node_farther
        logical, intent(in) :: lPressure
        logical :: lsumup, loverlap, lLeft
        integer :: dist_calc_num
        real :: dist_square
        integer :: next_cut_dim 

        ! if overlap, then return
        if (loverlap) return

        ! calculate the distance and include this point if it's within rcut
        ! only for the intermolecular part
        ! if .not. lsumup, account for chain number different from ichain
        ! if lsumup, only account for chain number greater than ichain
        if ((.not. lsumup .and. current_node%ichain .ne. ichain) .or. (lsumup .and. current_node%ichain .gt. ichain)) then
            dist_square = vector_dist_square(current_node%coord, ref_coord, rcut_square)
            dist_calc_num = dist_calc_num + 1
            if (dist_square .lt. rmin_square) then
                if (check_interaction(ichain, ibead, current_node%ichain, current_node%ibead)) then
                    loverlap = .true.
                    !write(io_output,*) "overlap in kdtree search: j, jj, dist_square", current_node%ichain, current_node%ibead, dist_square
                    !write(io_output,*) "overlap in kdtree search: coord", ref_coord, current_node%coord
                    return
                end if
            else if (dist_square .le. rcut_square) then
                search_output(1, current_dim) = dist_square
                search_output(2, current_dim) = current_node%ichain
                search_output(3, current_dim) = current_node%ibead

                if (lPressure) then
                    search_output(4, current_dim) = ref_coord(1) - current_node%coord(1)
                    search_output(5, current_dim) = ref_coord(2) - current_node%coord(2)
                    search_output(6, current_dim) = ref_coord(3) - current_node%coord(3)
                end if
                current_dim = current_dim + 1
            end if
        end if

        ! recursively search left or right in a more promising order
        next_cut_dim = cut_dim_next(cut_dim)
       
        if (l_left(ref_coord, current_node%coord, cut_dim)) then 
            node_closer => current_node%left_node
            node_farther => current_node%right_node
            lLeft = .true.
        else
            node_closer => current_node%right_node
            node_farther => current_node%left_node
            lLeft = .false.
        end if

        if (associated(node_closer)) then
            if (.not. node_closer%l_cube_updated) call set_cube(mol_tree, node_closer, lLeft)
            call range_search_in_tree(mol_tree, node_closer, ref_coord, next_cut_dim, ichain, ibead, lsumup, loverlap&
                    ,dist_calc_num, lPressure)
        end if

        ! search on the farther node if necessary
        if (associated(node_farther)) then
            if (.not. node_farther%l_cube_updated) call set_cube(mol_tree, node_farther, .not. lLeft)
            dist_square = dist_square_from_point_to_cube(mol_tree, ref_coord, node_farther, rcut_square)
            if (dist_square .le. rcut_square) then
                call range_search_in_tree(mol_tree, node_farther, ref_coord, next_cut_dim, ichain, ibead, lsumup, loverlap&
                        ,dist_calc_num, lPressure)
            end if 
        end if
    end subroutine range_search_in_tree

    ! Merge subroutine used in the merge sort
    subroutine merge(A,AINDEX,AA,AAA,NA,B,BINDEX,BB,BBB,NB,C,CINDEX,CC,CCC,NC)
        integer, intent(in) :: NA,NB,NC         ! Normal usage: NA+NB = NC
        real, intent(in out) :: A(NA), AA(NA), AAA(NA)
        integer, intent(in out) :: AINDEX(NA)        ! B overlays C(NA+1:NC)
        real, intent(in)     :: B(NB), BB(NB), BBB(NB)
        integer, intent(in out) :: BINDEX(NB)
        real, intent(in out) :: C(NC), CC(NC), CCC(NC)
        integer, intent(in out) :: CINDEX(NC)
 
        integer :: I,J,K
 
        I = 1; J = 1; K = 1;
        do while(I <= NA .and. J <= NB)
            if (A(I) < B(J) .or. (A(I).eq.B(J).and.AA(I).lt.BB(J)) .or. (A(I).eq.B(J).and.AA(I).eq.BB(J).and.AAA(I).lt.BBB(J))) then
                C(K) = A(I)
                CC(K) = AA(I)
                CCC(K) = AAA(I)
                CINDEX(K) = AINDEX(I)
                I = I+1
            else
                C(K) = B(J)
                CC(K) = BB(J)
                CCC(K) = BBB(J)
                CINDEX(K) = BINDEX(J)
                J = J+1
            endif
            K = K + 1
        enddo
        do while (I <= NA)
            C(K) = A(I)
            CC(K) = AA(I)
            CCC(K) = AAA(I)
            CINDEX(K) = AINDEX(I)
            I = I + 1
            K = K + 1
        enddo
        return
    end subroutine merge

    ! merge_sort subroutine
    recursive subroutine merge_sort(A,AINDEX,AA,AAA,N)
        integer, intent(in) :: N
        real, dimension(N), intent(in out) :: A,AA,AAA
        integer, dimension(N), intent(in out) :: AINDEX
        real, dimension((N+1)/2) :: T,TT,TTT
        integer, dimension((N+1)/2) :: TINDEX

        integer :: NA,NB,V
        real :: VA
 
        if (N < 2) return
        if (N == 2) then
            if (A(1).gt.A(2) .or. (A(1).eq.A(2).and.AA(1).gt.AA(2)) .or. (A(1).eq.A(2).and.AA(1).eq.AA(2).and.AAA(1).gt.AAA(2))) then
                VA = A(1)
                A(1) = A(2)
                A(2) = VA
                VA = AA(1)
                AA(1) = AA(2)
                AA(2) = VA
                VA = AAA(1)
                AAA(1) = AAA(2)
                AAA(2) = VA
                V = AINDEX(1)
                AINDEX(1) = AINDEX(2)
                AINDEX(2) = V
            end if
            return
        end if      
        NA=(N+1)/2
        NB=N-NA
 
        call merge_sort(A,AINDEX,AA,AAA,NA)
        call merge_sort(A(NA+1),AINDEX(NA+1),AA(NA+1),AAA(NA+1),NB)
 
        if (A(NA).gt.A(NA+1) .or. (A(NA).eq.A(NA+1).and.AA(NA).gt.AA(NA+1)) .or. &
            (A(NA).eq.A(NA+1).and.AA(NA).eq.AA(NA+1).and.AAA(NA).gt.AAA(NA+1))) then
            T(1:NA)=A(1:NA)
            TT(1:NA)=AA(1:NA)
            TTT(1:NA)=AAA(1:NA)
            TINDEX(1:NA)=AINDEX(1:NA)
            call merge(T,TINDEX,TT,TTT,NA,A(NA+1),AINDEX(NA+1),AA(NA+1),AAA(NA+1),NB,A,AINDEX,AA,AAA,N)
        endif
        return
    end subroutine merge_sort

    ! find the median of the coord_array
    ! sort the coord_array, chain_array and bead_array in ascending order of the coord_array
    ! coord_array: array consisting of coordinates
    ! N: the coord_array dimension
    ! median_index: returned value of the index of the median found
    ! index_array: returned array with the sorted index
    ! median_index_unsorted: returned avlue of the index of the median before sorting
    subroutine find_median(coord_array1, coord_array2, coord_array3, N, median_index, index_array, median_index_unsorted)
        integer, intent(in) :: N
        real, dimension(N), intent(in out) :: coord_array1, coord_array2, coord_array3
        integer, dimension(N) :: index_array
        integer, intent(out) :: median_index, median_index_unsorted
        integer :: i

        ! initialize the index_array
        do i = 1, N
            index_array(i) = i
        end do

        ! sort
        call merge_sort(coord_array1, index_array, coord_array2, coord_array3, N)
        
        ! find median
        if (mod(N, 2) == 0) then
            median_index = N / 2
        else
            median_index = (N + 1) / 2
        end if

        median_index_unsorted = index_array(median_index)

    end subroutine find_median

    ! Initialize the tree coordinates by inserting the median coordinates of the cutting dimension
    ! mol_tree: the pointer to the tree
    ! x_coord, y_coord, z_coord: N-dimensional array of x, y or z coordinates to be input
    ! chain_array, bead_array: N-dimensional array containing chain and bead number information, one-to-one correspondence to x_coord
    ! ix_array, iy_array, iz_array: N-dim array containing ix, iy and iz information
    ! Ndim: the dimension of the x_coord, y_coord and z_coord arrays
    ! N: the number of total beads
    ! cut_dim: the cutting dimension of the current node 
    recursive subroutine init_tree_coordinates(mol_tree, x_coord, y_coord, z_coord, chain_array, bead_array &
                            , ix_array, iy_array, iz_array, Ndim, N, cut_dim)
        type(tree), pointer :: mol_tree
        integer, intent(in) :: Ndim, N
        real, dimension(Ndim) :: x_coord, y_coord, z_coord
        real, dimension(N) :: coord_sort1, coord_sort2, coord_sort3
        real, allocatable :: x_coord_left(:), x_coord_right(:), y_coord_left(:), y_coord_right(:), z_coord_left(:), z_coord_right(:)
        integer, allocatable :: chain_array_left(:), chain_array_right(:), bead_array_left(:), bead_array_right(:)
        integer, allocatable :: ix_array_left(:), iy_array_left(:), iz_array_left(:)
        integer, allocatable :: ix_array_right(:), iy_array_right(:), iz_array_right(:)
        integer, dimension(Ndim) :: chain_array, bead_array, ix_array, iy_array, iz_array
        integer :: cut_dim, next_cut_dim, median_index, median_index_unsorted
        integer, dimension(N) :: index_array
        integer :: i, j, k, ichain, ibead, ix, iy, iz, nLeft, nRight
        real, dimension(3) :: coord

        if (cut_dim .eq. 1) then
            coord_sort1(1:N) = x_coord(1:N)
            coord_sort2(1:N) = y_coord(1:N)
            coord_sort3(1:N) = z_coord(1:N)
            next_cut_dim = 2
        else if (cut_dim .eq. 2) then
            coord_sort1(1:N) = y_coord(1:N)
            coord_sort2(1:N) = z_coord(1:N)
            coord_sort3(1:N) = x_coord(1:N)
            next_cut_dim = 3
        else
            coord_sort1(1:N) = z_coord(1:N)
            coord_sort2(1:N) = x_coord(1:N)
            coord_sort3(1:N) = y_coord(1:N)
            next_cut_dim = 1
        end if 

        ! find and insert the median node
        call find_median(coord_sort1, coord_sort2, coord_sort3, N, median_index, index_array, median_index_unsorted)
        coord(1) = x_coord(median_index_unsorted)
        coord(2) = y_coord(median_index_unsorted) 
        coord(3) = z_coord(median_index_unsorted)
        ichain = chain_array(median_index_unsorted)
        ibead = bead_array(median_index_unsorted)
        ix = ix_array(median_index_unsorted)
        iy = iy_array(median_index_unsorted)
        iz = iz_array(median_index_unsorted)
        mol_tree => insert_node(mol_tree, coord, ichain, ibead, ix, iy, iz)
        
        nLeft = median_index - 1
        nRight = N - median_index

        ! Split arrays for the recursive call
        if (nLeft .ge. 1) then
            allocate(x_coord_left(nLeft), y_coord_left(nLeft), z_coord_left(nLeft), chain_array_left(nLeft) &
                        ,bead_array_left(nLeft), ix_array_left(nLeft), iy_array_left(nLeft), iz_array_left(nLeft))
            do i = 1, nLeft
                k = index_array(i)
                x_coord_left(i) = x_coord(k)
                y_coord_left(i) = y_coord(k)
                z_coord_left(i) = z_coord(k)
                chain_array_left(i) = chain_array(k)
                bead_array_left(i) = bead_array(k)
                ix_array_left(i) = ix_array(k)
                iy_array_left(i) = iy_array(k)
                iz_array_left(i) = iz_array(k)
            end do
            
            call init_tree_coordinates(mol_tree, x_coord_left, y_coord_left, z_coord_left, &
                chain_array_left, bead_array_left, ix_array_left, iy_array_left, iz_array_left, nLeft, nLeft, next_cut_dim) 
        end if

        if (nRight .ge. 1) then
             allocate(x_coord_right(nRight), y_coord_right(nRight), z_coord_right(nRight), chain_array_right(nRight) &
                        , bead_array_right(nRight), ix_array_right(nRight), iy_array_right(nRight), iz_array_right(nRight))
            j = 1
            do i = median_index + 1, N
                k = index_array(i)
                x_coord_right(j) = x_coord(k)
                y_coord_right(j) = y_coord(k)
                z_coord_right(j) = z_coord(k)
                chain_array_right(j) = chain_array(k)
                bead_array_right(j) = bead_array(k)
                ix_array_right(j) = ix_array(k)
                iy_array_right(j) = iy_array(k)
                iz_array_right(j) = iz_array(k)
                j = j + 1
            end do
    
            call init_tree_coordinates(mol_tree, x_coord_right, y_coord_right, z_coord_right, &
                chain_array_right, bead_array_right, ix_array_right, iy_array_right, iz_array_right, nRight, nRight, next_cut_dim) 
        end if
         
    end subroutine init_tree_coordinates


    ! To construct/re-construct the kdtree
    !< ibox: for which simulation box
    !< iTree: if construct from scratch, iTree=ibox; if re-construct for volume move, iTree=nbox+1 or nbox+2
    !< lOutput: output the tree construction output if .true.
    subroutine construct_kdtree(ibox, iTree, lOutput)
        use sim_system
        use util_runtime,only:err_exit

        ! kd-tree variables
        logical :: lOutput, lAdd
        integer, intent(in) :: ibox, iTree
        type(interval), pointer :: cubeP(:)
        type(tree), pointer :: tree
        real, allocatable :: rxu_tot(:), ryu_tot(:), rzu_tot(:)
        integer, allocatable :: bead_array(:), chain_array(:), ix_array(:), iy_array(:), iz_array(:)
        integer :: i, ix, iy, iz, ichain, ibead, Ntot, Nactual
        integer, allocatable :: nx(:), ny(:), nz(:)
        integer :: imolty
        real :: xmin, xmax, ymin, ymax, zmin, zmax, xcoord, ycoord, zcoord, rbcut_plus_buffer

        if (.not. lkdtree_box(ibox)) call err_exit(__FILE__,__LINE__,'Cannot create kdtree for a box with lkdtree_box=F ',myid+1)

        Ntot = 0
        rbcut_plus_buffer = rcut(ibox) + kdtree_buffer_len(ibox)
        do ichain = 1, nchain
            if (nboxi(ichain) .eq. ibox) then
                imolty = moltyp(ichain)
                Ntot = Ntot + nunit(imolty)
            end if
        end do
        if (lpbc) Ntot = Ntot * 27

        if (lOutput) write(io_output, *) "Starting to construct the kd-tree"
        allocate(cubeP(3))
        do i = 1, 3
            cubeP(i)%lower = -2.5*boxlx(ibox)
            cubeP(i)%upper = 2.5*boxlx(ibox)
        end do

        if (associated(mol_tree(iTree)%tree)) call empty_tree(mol_tree(iTree)%tree)
        tree => init_tree(mol_tree(iTree)%tree, ibox, cubeP)

        if (lpbc) then
            allocate(nx(3), ny(3), nz(3))
            nx = [0, -1, 1]
            ny = [0, -1, 1]
            nz = [0, -1, 1]
        else
            allocate(nx(1), ny(1), nz(1))
            nx = [0]
            ny = [0]
            nz = [0]
        end if
    
        allocate(rxu_tot(Ntot), ryu_tot(Ntot),rzu_tot(Ntot), chain_array(Ntot), bead_array(Ntot) &
            ,ix_array(Ntot), iy_array(Ntot), iz_array(Ntot))

        i = 1
        xmin = 0.0E0_dp
        xmax = boxlx(ibox)
        ymin = 0.0E0_dp
        ymax = boxly(ibox)
        zmin = 0.0E0_dp
        zmax = boxlz(ibox)
        Nactual = 0

        do ix = 1, size(nx)
            do iy = 1, size(ny)
                do iz = 1, size(nz)
                    do ichain = 1, nchain
                        if (nboxi(ichain) .eq. ibox) then
                            imolty = moltyp(ichain)
                            do ibead = 1, nunit(imolty)
                                xcoord = rxu(ichain, ibead) + nx(ix) * boxlx(ibox)
                                ycoord = ryu(ichain, ibead) + ny(iy) * boxly(ibox)
                                zcoord = rzu(ichain, ibead) + nz(iz) * boxlz(ibox)

                                if ((ix .eq. 1) .and. (iy .eq. 1) .and. (iz .eq. 1)) then
                                    if (xcoord .gt. xmax) xmax = xcoord
                                    if (ycoord .gt. ymax) ymax = ycoord
                                    if (zcoord .gt. zmax) zmax = zcoord
                                    if (xcoord .lt. xmin) xmin = xcoord
                                    if (ycoord .lt. ymin) ymin = ycoord
                                    if (zcoord .lt. zmin) zmin = zcoord
                                        lAdd = .true.
                                else
                                    ! insert the coordinates only if it is within min-cutoff and max+cutoff
                                    if ((xcoord .gt. (xmin-rbcut_plus_buffer)) .and. (xcoord .lt. (xmax+rbcut_plus_buffer)) &
                                        .and. (ycoord .gt. (ymin-rbcut_plus_buffer)) .and. (ycoord .lt. (ymax+rbcut_plus_buffer)) &
                                        .and. (zcoord .gt. (zmin-rbcut_plus_buffer)) .and. (zcoord .lt. (zmax+rbcut_plus_buffer))) then
                                        lAdd = .true.
                                    else
                                        lAdd = .false.
                                    end if
                                end if

                                if (lAdd) then
                                    rxu_tot(i) = xcoord
                                    ryu_tot(i) = ycoord
                                    rzu_tot(i) = zcoord
                                    chain_array(i) = ichain
                                    bead_array(i) = ibead
                                    ix_array(i) = nx(ix)
                                    iy_array(i) = ny(iy)
                                    iz_array(i) = nz(iz)
                                    i = i + 1
                                    Nactual = Nactual + 1
                                end if
                            end do
                        end if
                    end do
                end do
            end do
        end do

        tree%bound(1)%lower = xmin
        tree%bound(1)%upper = xmax
        tree%bound(2)%lower = ymin
        tree%bound(2)%upper = ymax
        tree%bound(3)%lower = zmin
        tree%bound(3)%upper = zmax
        call init_tree_coordinates(tree, rxu_tot, ryu_tot, rzu_tot, chain_array, bead_array &
                , ix_array, iy_array, iz_array, Ntot, Nactual, 1)
        tree_height(ibox) = tree%height

        !write(io_output,*) 'Finished constructing the tree at ',time_date_str()
        if (lOutput) then
            write(io_output, *) "The number of nodes in the tree", tree%node_num
            write(io_output, *) "The height of the tree", tree%height
        end if

    end subroutine construct_kdtree
  
    ! allocate kdtree 
    subroutine allocate_kdtree()
        use sim_system
        integer :: ibox

        ! allocate allocatable arrays
        if (allocated(mol_tree)) deallocate(mol_tree, tree_height, lkdtree_box, kdtree_buffer_len)
        allocate(mol_tree(nbox+2), tree_height(nbox+2)) !< nbox+2 because the extra 2 will be used to temporarily store tree in volume moves
        allocate(lkdtree_box(nbox),kdtree_buffer_len(nbox))

        ! set default value of lkdtree_box to be false
        do ibox = 1, nbox
            lkdtree_box(ibox) = .false.
        end do

        ! reset tree pointers
        do ibox = 1, nbox+2
            if (associated(mol_tree(ibox)%tree)) nullify(mol_tree(ibox)%tree)
        end do
 
    end subroutine allocate_kdtree


    ! update the kdtree for the volume move
    ! point mol_tree(ibox)%tree to the mol_tree(nbox+1)%tree
    subroutine update_box_kdtree(ibox)
        use sim_system

        integer, intent(in) :: ibox
        integer :: iTree !< the index of the updated mol_tree
        integer :: i

        ! find iTree
        iTree = 0
        do i = nbox+1, nbox+2
            if (associated(mol_tree(i)%tree)) then
                if (mol_tree(i)%tree%box .eq. ibox) then
                    iTree = i
                    exit
                end if
            end if
        end do 
        if (iTree .eq. 0) call err_exit(__FILE__,__LINE__,'Error in update_box_kdtree: iTree not found',myid)

        ! Empty the old tree 
        call empty_node(mol_tree(ibox)%tree%tree_root)
 
        ! Update the new tree
        mol_tree(ibox)%tree%tree_root => mol_tree(iTree)%tree%tree_root
        mol_tree(ibox)%tree%height = mol_tree(iTree)%tree%height
        mol_tree(ibox)%tree%node_num = mol_tree(iTree)%tree%node_num 
        mol_tree(ibox)%tree%cube = mol_tree(iTree)%tree%cube 
        mol_tree(ibox)%tree%bound = mol_tree(iTree)%tree%bound
        mol_tree(ibox)%tree%box = mol_tree(iTree)%tree%box
        deallocate(mol_tree(iTree)%tree) 
        
    end subroutine update_box_kdtree

    ! read the kdtree related variables
    subroutine read_kdtree(io_input)
        use var_type,only:default_path_length,default_string_length
        use util_string,only:uppercase
        use util_files,only:get_iounit,readLine
        use util_mp,only:mp_bcast
        use sim_system
        integer,intent(in)::io_input
        character(LEN=default_string_length)::line_in
        integer :: jerr

        ! Check compatability of other settings with KDTREE
        if (lijall .and. lkdtree) call err_exit(__FILE__,__LINE__,'KDTREE cannot be used when lijall is true',myid+1)        
        if (lchgall .and. lkdtree) call err_exit(__FILE__,__LINE__,'KDTREE cannot be used when lchgall is true',myid+1)        
        if (lneigh .and. lkdtree) call err_exit(__FILE__,__LINE__,'KDTREE cannot be used when neighbor list is used',myid+1)        

        ! Look for section KDTREE
        if ((myid .eq. rootid) .and. (lkdtree)) then
            REWIND(io_input)
            CYCLE_READ_KDTREE:DO
                call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Section KDTREE not found',jerr)
                if (UPPERCASE(line_in(1:10)).eq.'KDTREE') then
                    exit cycle_read_kdtree
                end if
            END DO CYCLE_READ_KDTREE
        
            call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
            if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section KDTREE',jerr)
            if (UPPERCASE(line_in(1:14)).eq.'END KDTREE') call err_exit(__FILE__,__LINE__&
                        ,'Section KDTREE not complete!',myid+1)
            read(line_in,*) lkdtree_box(1:nbox)
            
            call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
            if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section KDTREE',jerr)
            if (UPPERCASE(line_in(1:14)).eq.'END KDTREE') call err_exit(__FILE__,__LINE__&
                        ,'Section KDTREE not complete!',myid+1)
            read(line_in,*) kdtree_buffer_len(1:nbox)

        end if

        call mp_bcast(lkdtree_box,nbox,rootid,groupid)
        call mp_bcast(kdtree_buffer_len,nbox,rootid,groupid)

    end subroutine read_kdtree

    function check_interaction(i, ii, j, jj) result(l_interaction)
        use sim_system,only:moltyp,ntype,lij,lqchg

        integer, intent(in) :: i, ii, j, jj
        logical :: l_interaction
        integer :: imolty, jmolty, ntii, ntjj
        logical :: lij_i, lij_j, qq_i, qq_j

        imolty = moltyp(i)
        jmolty = moltyp(j)
        ntii = ntype(imolty, ii)
        ntjj = ntype(jmolty, jj)
        lij_i = lij(ntii)
        lij_j = lij(ntjj)
        qq_i = lqchg(ntii)
        qq_j = lqchg(ntjj)

        if ((.not. (lij_i .and. lij_j)) .and. (.not. (qq_i .and. qq_j))) then
            l_interaction = .false.
        else
            l_interaction = .true.
        end if

        return
    end function check_interaction
END MODULE util_kdtree
