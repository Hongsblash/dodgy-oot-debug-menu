.include "macro.inc"

# assembler directives
.set noat      # allow manual use of $at
.set noreorder # don't insert nops after branches
.set gp=64     # allow use of 64-bit general purpose registers

.section .data

.balign 16

glabel D_8012DBA0
    .incbin "baserom.z64", 0xBA4D40, 0x4

glabel D_8012DBA4
    .incbin "baserom.z64", 0xBA4D44, 0xC
