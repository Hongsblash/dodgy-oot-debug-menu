.include "macro.inc"

 # assembler directives
 .set noat      # allow manual use of $at
 .set noreorder # don't insert nops after branches
 .set gp=64     # allow use of 64-bit general purposee registers

.section .data

glabel Boss_Tw_InitVars
 .word 0x00DC0900, 0x00000035, 0x00D30000, 0x000006B4
.word BossTw_Init
.word BossTw_Destroy
.word BossTw_Update
.word BossTw_Draw
glabel D_8094A7D0
 .word 0x00000000, 0x00000000, 0x447A0000
glabel D_8094A7DC
 .word 0x00000000, 0x00000000, 0x00000000
glabel D_8094A7E8
 .word 0x0A390909, 0x10010000, 0x00000000, 0xFFCFFFFF, 0x00300000, 0x00100000, 0x00000000, 0x01010100, 0x00190023, 0xFFEF0000, 0x00000000
glabel D_8094A814
 .word 0x03110909, 0x10010000, 0x00000000, 0xFFCFFFFF, 0x00200000, 0xFFCDFFFE, 0x00000000, 0x01010100, 0x002D0078, 0xFFE20000, 0x00000000
glabel D_8094A840
 .word 0x03110939, 0x10010000, 0x00000000, 0xFFCFFFFF, 0x00200000, 0xFFCDFFFE, 0x00000000, 0x01050100, 0x002D0078, 0xFFE20000, 0x00000000
glabel D_8094A86C
 .word 0x44110000, 0x43BE0000, 0x00000000, 0x00000000, 0x43BE0000, 0x44110000, 0xC4110000, 0x43BE0000, 0x00000000, 0x00000000, 0x43BE0000, 0xC4110000
glabel D_8094A89C
 .word 0x00000000
glabel D_8094A8A0
 .word 0x801F0005, 0xB86C0000, 0x304C0000
glabel D_8094A8AC
 .word 0x44160000, 0x43C80000, 0x00000000, 0x00000000, 0x43C80000, 0x44160000, 0xC4160000, 0x43C80000, 0x00000000, 0x00000000, 0x43C80000, 0xC4160000
glabel D_8094A8DC
 .word 0x00000000, 0x00000000, 0x00000000
glabel D_8094A8E8
 .word 0x00000000, 0x00000000, 0x00000000
glabel D_8094A8F4
 .word 0x00000000, 0x00000000, 0x00000000
glabel D_8094A900
 .word 0x00000001, 0x00020002, 0x00010000
glabel D_8094A90C
 .word 0x00000001, 0x00020002, 0x00020002, 0x00020002, 0x00010000
glabel D_8094A920
 .word 0x00000000, 0x00000000, 0x00000000
glabel D_8094A92C
 .word 0x00000000, 0x00000000, 0x00000000
glabel D_8094A938
 .word 0x00000000, 0x00000000, 0x00000000
glabel D_8094A944
 .word 0x00000000, 0x00000000, 0x00000000
glabel D_8094A950
 .word 0x00000000, 0x44FA0000, 0xC4FA0000
glabel D_8094A95C
 .word 0x00000000, 0x00000000, 0xC61C4000
glabel D_8094A968
 .word 0x00000000, 0x00000000, 0xC5FA0000
glabel D_8094A974
 .word 0x00000000, 0x00000000, 0xC60CA000
glabel D_8094A980
 .word 0x00000000, 0x00000000, 0xC62BE000
glabel D_8094A98C
 .word 0x00000000, 0x00000000, 0xC63B8000
glabel D_8094A998
 .word 0x0600A438, 0x0600B238, 0x0600B638
glabel D_8094A9A4
 .word 0x00000000, 0x43480000, 0x44FA0000
glabel D_8094A9B0
 .word 0x0602A9B0, 0x0602A070, 0x0602A470
glabel D_8094A9BC
 .word 0x00000000, 0x00000000, 0x00000000
glabel D_8094A9C8
 .word 0x00000000, 0x44FA0000, 0xC4FA0000
glabel D_8094A9D4
 .word 0x464B2000, 0x00000000, 0x00000000
glabel D_8094A9E0
 .word 0x464B2000, 0x00000000, 0x00000000
glabel D_8094A9EC
 .word 0x00000000, 0x43480000, 0x44FA0000
glabel D_8094A9F8
 .word 0x00000000, 0x00000000, 0x00000000
glabel D_8094AA04
 .word 0x00000000, 0x00000000, 0x00000000
glabel D_8094AA10
 .word 0x00000000, 0x00000000, 0x00000000
glabel D_8094AA1C
 .word 0x00000000, 0x00000000, 0x00000000
glabel D_8094AA28
 .word 0xFF8000FF, 0x0000FFFF, 0x00FF0000, 0x646464FF, 0xFFFF9696, 0x96FFFFFF
glabel D_8094AA40
 .word 0x00000000, 0x00000000, 0x00000000
glabel D_8094AA4C
 .word 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000

