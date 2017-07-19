module nesting_support
    use params_common_sn
    use communication_helper_mpi
    implicit none
contains
    subroutine calcSubgridCoords(rank,i_s,j_s)
        integer, intent(In) :: rank
        integer, intent(Out) :: i_s,j_s
        integer, dimension(2) :: ij_s
        integer :: ierr
        call MPI_Cart_coords(cartTopComm, rank, 2, ij_s, ierr)
        i_s = ij_s(1)
        j_s = ij_s(2)
    end subroutine calcSubgridCoords

    subroutine currentSubgridCoords(i_s,j_s)
        integer, intent(Out) :: i_s,j_s
        integer, dimension(2) :: ij_s
        integer :: ierror, rank
        call MPI_COMM_Rank(communicator, rank, ierror)
        call checkMPIError()
        call MPI_Cart_coords(cartTopComm, rank, 2, ij_s, ierror)
        call checkMPIError()
        i_s = ij_s(1)
        j_s = ij_s(2)
    end subroutine currentSubgridCoords

    logical function inNestedGridByRank(rank) result(in_grid)
            integer, intent(In) :: rank
            integer :: i_s, j_s
            call calcSubgridCoords(rank,i_s,j_s)
            in_grid = i_s >= i_s_nest_start .and. i_s <= i_s_nest_end .and. j_s >= j_s_nest_start .and. j_s <= j_s_nest_end
    end function inNestedGridByRank

    logical function inNestedGrid() result(in_grid)
            integer :: rank ,ierror
            call MPI_COMM_Rank(communicator, rank, ierror)
            call checkMPIError()
            in_grid = inNestedGridByRank(rank)
    end function inNestedGrid

end module nesting_support
