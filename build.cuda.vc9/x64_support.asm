
    BITS 64

;   uint64 mul_mod(uint64 x, uint64 y, uint64 m)
;   return (x * y) % m
     
    global mul_mod_64

mul_mod_64:
    mov     rax, rcx
    mul     rdx
    div     r8
    mov     rax, rdx
    ret
