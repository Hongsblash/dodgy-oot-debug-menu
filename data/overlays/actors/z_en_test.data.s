.include "macro.inc"

 # assembler directives
 .set noat      # allow manual use of $at
 .set noreorder # don't insert nops after branches
 .set gp=64     # allow use of 64-bit general purposee registers

.section .data

glabel D_80864510
 .word 0x00000000, 0x01010101, 0x01010101, 0x01010101, 0x01010101, 0x01010101, 0x01010101, 0x01010101, 0x01010101, 0x01010100, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000
glabel En_Test_InitVars
 .word 0x00020500, 0x00000015, 0x00320000, 0x00000928
.word EnTest_Init
.word EnTest_Destroy
.word EnTest_Update
.word EnTest_Draw
glabel D_80864570
 .word 0x05000939, 0x10010000, 0x00000000, 0x00000000, 0x00000000, 0xFFCFFFFF, 0x00000000, 0x00010100, 0x00190041, 0x00000000, 0x00000000
glabel D_8086459C
 .word 0x09000D00, 0x00010000, 0x00000000, 0x00000000, 0x00000000, 0xFFC1FFFF, 0x00000000, 0x00010000, 0x00140046, 0xFFCE0000, 0x00000000
glabel D_808645C8
 .word 0x0A110000, 0x00030000, 0x00000000, 0xFFCFFFFF, 0x00100000, 0x00000000, 0x00000000, 0x81000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000
glabel D_80864618
 .word 0x1002D102, 0x10020210, 0x01020402, 0xF4E20202, 0x0260F3E0, 0x00000104, 0x02020804, 0x00000400
glabel D_80864638
 .word 0x8917001B, 0xB04C01F4, 0xC850000F, 0xB0540000, 0x386CFA24, 0x44898000, 0xC42F0000, 0x00000000
glabel D_80864658
 .word 0x43960000, 0x00000000, 0x00000000
glabel D_80864664
 .word 0x45548000, 0x00000000, 0x00000000
glabel D_80864670
 .word 0x00000000, 0x00000000, 0x00000000
glabel D_8086467C
 .word 0x45DAC000, 0x447A0000, 0x00000000
glabel D_80864688
 .word 0x453B8000, 0xC4FA0000, 0xC47A0000
glabel D_80864694
 .word 0x453B8000, 0xC4FA0000, 0x447A0000
glabel D_808646A0
 .word 0xC4A28000, 0x44898000, 0x00000000, 0xC53B8000, 0x44ED8000, 0x44480000, 0xC53B8000, 0xC4898000, 0x44480000, 0x44ED8000, 0x44ED8000, 0x44480000, 0xC53B8000, 0xC4898000, 0x44480000, 0x44ED8000, 0xC4898000, 0x44480000, 0x44ED8000, 0x44ED8000, 0x44480000, 0x00000000, 0x00000000, 0x00000000

