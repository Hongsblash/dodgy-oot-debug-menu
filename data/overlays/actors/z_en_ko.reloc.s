.include "macro.inc"

 # assembler directives
 .set noat      # allow manual use of $at
 .set noreorder # don't insert nops after branches
 .set gp=64     # allow use of 64-bit general purposee registers

.section .rodata
glabel D_80A9A9E0

.incbin "baserom/ovl_En_Ko", 0x3C30, 0x000005C0
