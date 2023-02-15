subroutine vert2cell(nEdgesOnCell, verticesOnCell, maxEdges, nVertLevels, nCells, nVertices, fieldv, fieldc)
implicit none

integer, intent(in) :: maxEdges, nVertLevels, nCells, nVertices
!f2py required :: maxEdges, nVertLevels, nCells, nVertices
integer, intent(in) :: nEdgesOnCell(nCells)
integer, intent(in) :: verticesOnCell(maxEdges,nCells)
real*8, intent(in)  :: fieldv(nVertLevels,nVertices)
real*8, intent(out) :: fieldc(nVertLevels,nCells)

integer i,j,k
real*8 factor

do k=1,nVertLevels
    do i=1,nCells
        factor = 1./nEdgesOnCell(i)
        fieldc(k,i) = 0.
        do j=1,nEdgesOnCell(i)
            fieldc(k,i)=fieldc(k,i) + fieldv(k,verticesOnCell(j,i))
        end do
            fieldc(k,i) = factor*fieldc(k,i)
    end do
end do
return
end subroutine vert2cell
