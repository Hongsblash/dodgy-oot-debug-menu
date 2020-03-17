.include "macro.inc"

 # assembler directives
 .set noat      # allow manual use of $at
 .set noreorder # don't insert nops after branches
 .set gp=64     # allow use of 64-bit general purposee registers

.section .data

glabel D_80AB85E0
 .word 0x00000000
glabel En_Niw_InitVars
 .word 0x00190600, 0x00800010, 0x00130000, 0x000007B8
.word EnNiw_Init
.word EnNiw_Destroy
.word EnNiw_Update
.word EnNiw_Draw
glabel D_80AB8604
 .word 0x459C4000, 0xC59C4000
glabel D_80AB860C
 .word 0x459C4000
glabel D_80AB8610
 .word 0x453B8000, 0x457A0000
glabel D_80AB8618
 .word 0xC4D42000, 0x42A00000, 0x44598000, 0x42640000, 0x43A00000, 0xC4284000, 0x44470000, 0x42A00000, 0x44CCE000, 0x44B12000, 0x43E88000, 0x43290000, 0xC2700000, 0x00000000, 0xC2380000, 0xC3770000, 0x42A00000, 0x44558000, 0x4486E000, 0x42A00000, 0xC23C0000
glabel D_80AB866C
 .word 0x02000400, 0x08001000, 0x20004000, 0x80000000
glabel D_80AB867C
 .word 0x00000000
glabel D_80AB8680
 .word 0x00000000
glabel D_80AB8684
 .word 0x05000901, 0x20010000, 0x00000000, 0x00000000, 0x00000000, 0xFFCFFFFF, 0x00000000, 0x00010100, 0x000F0019, 0x00040000, 0x00000000
glabel D_80AB86B0
 .word 0x0A000039, 0x20010000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000100, 0x000F0019, 0x00040000, 0x00000000
glabel D_80AB86DC
 .word 0x801F0006, 0xB86CF830, 0x304C0000
glabel D_80AB86E8
 .word 0x00000000, 0x00000000, 0x00000000
glabel D_80AB86F4
 .word 0x00000000, 0x00000000, 0x00000000
glabel D_80AB8700
 .word 0x00000000, 0x00000000, 0x00000000
glabel D_80AB870C
 .word 0x3E19999A, 0x3E19999A, 0x3E19999A, 0x00000000, 0x00000000

