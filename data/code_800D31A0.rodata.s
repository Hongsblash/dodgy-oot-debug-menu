.include "macro.inc"

# assembler directives
.set noat      # allow manual use of $at
.set noreorder # don't insert nops after branches
.set gp=64     # allow use of 64-bit general purposee registers

.section .rodata

.align 4

glabel D_80146670
    .asciz "\x1b[31m\n**** Freeze!! ****\n\x1b[m"
    .balign 4
