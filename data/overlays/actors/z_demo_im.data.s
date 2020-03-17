.include "macro.inc"

 # assembler directives
 .set noat      # allow manual use of $at
 .set noreorder # don't insert nops after branches
 .set gp=64     # allow use of 64-bit general purposee registers

.section .data

glabel D_80987830
 .word 0x06007210, 0x06007D50, 0x06008150
glabel D_8098783C
 .word 0x00000000
glabel D_80987840
 .word 0x00000009, 0x01000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000100, 0x00190050, 0x00000000, 0x00000000
glabel D_8098786C
 .word 0x00000020, 0x00000BB8, 0x00000020, 0x00000001, 0x00010000, 0x0BB80000, 0x00000000, 0x00000000, 0xFFFFFFFC, 0x00000002, 0x00000000, 0xFFFFFFFC, 0x00000002, 0x00000000, 0x00000000, 0x00000000, 0x0000001F, 0x00000005, 0x00010000, 0x02B90000, 0x00000000, 0x00000000, 0x000000D8, 0x00000000, 0x00000000, 0x000000D8, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x000202B9, 0x02BA0000, 0x00000000, 0x00000000, 0x000000D8, 0x00000000, 0x00000000, 0x000000D8, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x000402BA, 0x03000000, 0x00000000, 0x00000000, 0x000000D8, 0x00000000, 0x00000000, 0x000000D8, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00020300, 0x03310000, 0x00000000, 0x00000000, 0x000000D8, 0x00000000, 0x00000000, 0x00000052, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00030331, 0x0A6A0000, 0x00000000, 0x00000000, 0x00000052, 0x00000000, 0x00000000, 0x00000052, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x0000000A, 0x00000003, 0x000D0000, 0x012C0000, 0x00000000, 0x00000000, 0x00000006, 0x00000000, 0x00000000, 0x00000006, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x0005012C, 0x02950000, 0xEAAA0000, 0x00000000, 0x00000006, 0x00000000, 0x00000000, 0x00000006, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00130295, 0x078E0000, 0x6AAA0000, 0x00000000, 0x00000006, 0x00000000, 0x00000000, 0x00000006, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x0000002C, 0x00000003, 0x00010000, 0x00910000, 0x00000000, 0xFFFFFF9F, 0x00000006, 0x000000A9, 0xFFFFFF9F, 0x00000006, 0x000000A9, 0x00000000, 0x00000000, 0x00000000, 0x00020091, 0x02670000, 0x00000000, 0xFFFFFF9F, 0x00000006, 0x000000A9, 0xFFFFFF9F, 0x00000006, 0x000000A9, 0x00000000, 0x00000000, 0x00000000, 0x00030267, 0x07720000, 0x00000000, 0xFFFFFF9F, 0x00000006, 0x000000A9, 0xFFFFFF9F, 0x00000006, 0x000000A9, 0x00000000, 0x00000000, 0x00000000, 0x00000031, 0x00000001, 0x00010000, 0x0BB80000, 0x00000000, 0xFFFFFFEA, 0x00000000, 0xFFFFFFC9, 0xFFFFFFEA, 0x00000000, 0xFFFFFFC9, 0x00000000, 0x00000000, 0x00000000, 0x00000004, 0x00000002, 0x00020000, 0x000A0000, 0x00000000, 0xFFFFFFFE, 0x00000000, 0x0000000D, 0xFFFFFFFE, 0x00000000, 0x0000000D, 0x00000000, 0x00000000, 0x00000000, 0x0002000A, 0x0BB80000, 0x00000000, 0xFFFFFFFE, 0x00000000, 0x0000000D, 0xFFFFFFFE, 0x00000000, 0x0000000D, 0x00000000, 0x00000000, 0x00000000, 0x0000002D, 0x00000001, 0x000502B6, 0x02D402D4, 0x0000002D, 0x00000001, 0x000103C0, 0x03DE03DE, 0x0000002D, 0x00000001, 0x000102AB, 0x02B402B4, 0x0000003E, 0x00000002, 0x00010000, 0x000A0000, 0x00000000, 0x00000040, 0x00000050, 0x00000082, 0x00000040, 0x00000050, 0x00000082, 0x00000000, 0x00000000, 0x00000000, 0x0004000A, 0x0BB80000, 0x00000000, 0x00000040, 0x00000050, 0x00000082, 0x00000040, 0x00000050, 0x00000082, 0x00000000, 0x00000000, 0x00000000, 0x00000056, 0x00000001, 0x00440302, 0x03030000, 0x00000000, 0x00000000, 0xFFFFFFC9, 0x0000005C, 0x00000000, 0xFFFFFFC9, 0x0000005C, 0x00000000, 0x00000000, 0x00000000, 0x00000013, 0x0000000E, 0xFFFF0000, 0x0154FFFF, 0xFFFFFFFF, 0x50220154, 0x01610000, 0x00000000, 0xFFFF0161, 0x0176FFFF, 0xFFFFFFFF, 0x50250176, 0x01940000, 0x00000000, 0xFFFF0194, 0x01A8FFFF, 0xFFFFFFFF, 0x502B01A8, 0x01DA0000, 0x00000000, 0xFFFF01DA, 0x01EEFFFF, 0xFFFFFFFF, 0x502C01EE, 0x021F0000, 0x00000000, 0xFFFF021F, 0x0234FFFF, 0xFFFFFFFF, 0x50260234, 0x02650000, 0x00000000, 0xFFFF0265, 0x03BBFFFF, 0xFFFFFFFF, 0x004103BB, 0x03BF0000, 0x00000000, 0xFFFF03BF, 0x03FCFFFF, 0xFFFFFFFF, 0x502303FC, 0x04050000, 0x00000000, 0x000003E8, 0x00000001, 0x00610424, 0x04250425, 0x0000007C, 0x00000001, 0x000402A1, 0x02D30000, 0x00000000, 0x00000000, 0xFFFFFFC1, 0x00000058, 0x00000000, 0xFFFFFFC1, 0x00000058, 0x00000000, 0x00000000, 0x00000000, 0x00000001, 0x00010000, 0x01550000, 0x00000000, 0x42726667, 0x002100E1, 0xFFC620BA, 0x00000000, 0x42726667, 0x002100E1, 0xFFC6D5E0, 0x00000000, 0x42726667, 0x002100E1, 0xFFC60950, 0x00000000, 0x42726667, 0x0021006A, 0xFFC67C50, 0x00000000, 0x42726667, 0x00210017, 0xFFC60000, 0x00000000, 0x42726667, 0x000B000A, 0xFFEEFFFF, 0x00000000, 0x42726667, 0x000B000A, 0xFFEE0000, 0x00000000, 0x42726667, 0x000B000A, 0xFFEEE6A0, 0xFF000000, 0x42726667, 0x000B000A, 0xFFEE7C53, 0x00000001, 0x00010107, 0x01F80000, 0x00000000, 0x41B50402, 0xFFCF000D, 0x009E20BA, 0x00000000, 0x41B50402, 0xFFCF000D, 0x009ED5E0, 0x00000000, 0x41B50402, 0xFFCF000D, 0x009E0950, 0x00000000, 0x41B50402, 0xFFEB0016, 0x00967C50, 0x00000000, 0x4204E872, 0xFFEB0016, 0x00960000, 0x00000000, 0x4204E872, 0xFFEB0016, 0x0096FFFF, 0x00000000, 0x4204E872, 0xFFEB0016, 0x00960000, 0x00000000, 0x4204E872, 0xFFEB0016, 0x0096E6A0, 0xFF000000, 0x4204E872, 0xFFEB0016, 0x00967C53, 0x00000001, 0x0001016B, 0x03380000, 0x00000000, 0x42726667, 0xFFBD0008, 0x007520BA, 0x00000000, 0x42726667, 0xFFBD0008, 0x0075D5E0, 0x00000000, 0x42726667, 0xFFBD0008, 0x00750950, 0x00000000, 0x42726667, 0xFFBD0008, 0x00757C50, 0x00000000, 0x42726667, 0xFFBD0008, 0x00750000, 0x00000000, 0x42726667, 0xFFBD0008, 0x0075FFFF, 0x00000000, 0x428D998E, 0xFFBD0008, 0x00750000, 0x00000000, 0x428D998E, 0xFFBD0008, 0x0075E6A0, 0x00000000, 0x428D998E, 0xFFBD0008, 0x00757C53, 0x00000000, 0x428D998E, 0xFFBD0008, 0x00750000, 0xFF000000, 0x428D998E, 0xFFBD0008, 0x00750000, 0x00000005, 0x0001019D, 0x05E00000, 0x00000000, 0x419FFFB1, 0xFFCB001C, 0x002D20BA, 0x00000000, 0x419FFFB1, 0xFFCB001C, 0x002DD5E0, 0x00000000, 0x419FFFB1, 0xFFCB001C, 0x002D0950, 0x00000000, 0x419FFFB1, 0xFFCB001C, 0x002D7C50, 0xFF000000, 0x419FFFB1, 0xFFCB001C, 0x002D0000, 0x00000001, 0x000101E3, 0x06940000, 0x00000000, 0x4289D68F, 0xFFBD0008, 0x007520BA, 0x00000000, 0x4289D68F, 0xFFBD0008, 0x0075D5E0, 0x00000000, 0x4289D68F, 0xFFBD0008, 0x00750950, 0x00000000, 0x4289D68F, 0xFFB9000F, 0x007C7C50, 0x00000000, 0x4289D68F, 0xFFB7001D, 0x00830000, 0x00000000, 0x4289D68F, 0xFFB20028, 0x008CFFFF, 0x00000000, 0x4289D68F, 0xFFB20028, 0x008C0000, 0x00000000, 0x4289D68F, 0xFFB20028, 0x008CE6A0, 0xFF000000, 0x4289D68F, 0xFFB20028, 0x008C7C53, 0x00000001, 0x00010229, 0x066C0000, 0x00000000, 0x42366658, 0xFFD90022, 0x00C920BA, 0x00000000, 0x42366658, 0xFFD90022, 0x00C9D5E0, 0x00000000, 0x42366658, 0xFFD90022, 0x00C90950, 0x00000000, 0x42366658, 0xFFD90022, 0x00C97C50, 0xFF000000, 0x42366658, 0xFFD90022, 0x00C90000, 0x00000001, 0x0001026F, 0x03330000, 0x00000000, 0x428D3328, 0x0009000D, 0xFFEF20BA, 0x00000000, 0x428D3328, 0x0009000D, 0xFFEFD5E0, 0x00000000, 0x428D3328, 0x0009000D, 0xFFEF0950, 0x00000000, 0x428D3328, 0x00090047, 0xFFEF7C50, 0x00000000, 0x42FDFF84, 0x00090181, 0xFFEF0000, 0x00000000, 0x42FDFF84, 0x00090181, 0xFFEFFFFF, 0x00000000, 0x42FDFF84, 0x00090181, 0xFFEF0000, 0x00000000, 0x42FDFF84, 0x00090181, 0xFFEFE6A0, 0xFF000000, 0x42FDFF84, 0x00090181, 0xFFEF7C53, 0x00000001, 0x000102B5, 0x040B0000, 0x00000000, 0x42700000, 0x000D0356, 0x000220BA, 0x00000000, 0x42700000, 0x00090355, 0x0005D5E0, 0x00000000, 0x42700000, 0xFFFD0355, 0x00050950, 0x00000000, 0x42700000, 0xFFF70355, 0xFFFA7C50, 0x00000000, 0x42700000, 0xFFFE0354, 0xFFEF0000, 0x00000000, 0x42700000, 0x00090354, 0xFFEFFFFF, 0x00000000, 0x42700000, 0x00100354, 0xFFFA0000, 0x00000000, 0x42700000, 0x00090354, 0x0005E6A0, 0xFF000000, 0x42700000, 0xFFFD0353, 0x00057C53, 0x00000005, 0x00010301, 0x079E0000, 0x00000000, 0x4289332C, 0x00000021, 0xFFE520BA, 0x00000000, 0x4289332C, 0x00000021, 0xFFE5D5E0, 0x00000000, 0x4289332C, 0x00000044, 0xFFE60950, 0x00000000, 0x4289332C, 0x00000067, 0xFFE67C50, 0x00000000, 0x4289332C, 0x00000067, 0xFFE60000, 0x00000000, 0x4289332C, 0x00000067, 0xFFE6FFFF, 0x00000000, 0x4289332C, 0x00000067, 0xFFE60000, 0xFF000000, 0x4289332C, 0x00000067, 0xFFE6E6A0, 0x00000002, 0x00010000, 0x01720000, 0x0000001E, 0x42726667, 0xFFEB0015, 0x002A20BA, 0x00000032, 0x42726667, 0xFFEB0015, 0x002AD5E0, 0x00000032, 0x42726667, 0xFFEB0015, 0x002A0950, 0x00000032, 0x42726667, 0xFFB20056, 0x00907C50, 0x00000032, 0x42726667, 0xFFB00015, 0x008E0000, 0x00000032, 0x42726667, 0xFF9D001F, 0x00B1FFFF, 0x0000001E, 0x42726667, 0xFF9D001F, 0x00B10000, 0x0000001E, 0x42726667, 0xFF9D001F, 0x00B1E6A0, 0xFF00001E, 0x42726667, 0xFF9D001F, 0x00B17C53, 0x00000002, 0x00010107, 0x02150000, 0x0000001E, 0x41B50402, 0xFECB000B, 0x00E520BA, 0x0000001E, 0x41B50402, 0xFECB000B, 0x00E5D5E0, 0x0000001E, 0x41B50402, 0xFECB000B, 0x00E50950, 0x0000001E, 0x4204E872, 0xFEFE0085, 0x00BF7C50, 0x0000001E, 0x4204E872, 0xFEFE0085, 0x00BF0000, 0x0000001E, 0x4204E872, 0xFEFE0085, 0x00BFFFFF, 0x0000001E, 0x4204E872, 0xFEFE0085, 0x00BF0000, 0x0000001E, 0x4204E872, 0xFEFE0085, 0x00BFE6A0, 0xFF00001E, 0x4204E872, 0xFEFE0085, 0x00BF7C53, 0x00000002, 0x0001016B, 0x03550000, 0x00000032, 0x42899992, 0x0044006D, 0xFF9920BA, 0x00000032, 0x4289FFF8, 0x0044006D, 0xFF99D5E0, 0x00000032, 0x42899992, 0x0043006D, 0xFF990950, 0x00000032, 0x4289332C, 0x00BD006D, 0x005C7C50, 0x00000032, 0x4289332C, 0x0039006C, 0x01550000, 0x00000032, 0x428D998E, 0xFF46006B, 0x0155FFFF, 0x00000032, 0x428D998E, 0xFF46006B, 0x01550000, 0x00000032, 0x428D998E, 0xFF46006B, 0x0155E6A0, 0x0000001E, 0x428D998E, 0xFF46006B, 0x01557C53, 0x0000001E, 0x428D998E, 0xFF46006B, 0x01550000, 0xFF00001E, 0x428D998E, 0xFF46006B, 0x01550000, 0x00000006, 0x0001019D, 0x05FD0000, 0x0000001E, 0x419FFFB1, 0x00AB007A, 0xFF9620BA, 0x0000001E, 0x419FFFB1, 0x00AB007A, 0xFF96D5E0, 0x000003E8, 0x419FFFB1, 0x00AB007A, 0xFF960950, 0x0000001E, 0x419FFFB1, 0x00AA007A, 0xFF967C50, 0xFF00001E, 0x419FFFB1, 0x00AA007A, 0xFF960000, 0x00000002, 0x000101E3, 0x06B10000, 0x00010014, 0x4289D68F, 0xFF410084, 0x014720BA, 0x00000014, 0x4289D68F, 0xFF410084, 0x0147D5E0, 0x00FF001E, 0x4289D68F, 0xFF420083, 0x01460950, 0x00000028, 0x4289D68F, 0xFF43009A, 0x01447C50, 0x0001001E, 0x4289D68F, 0xFF3E00A2, 0x014B0000, 0x0000001E, 0x4289D68F, 0xFF3900AC, 0x0153FFFF, 0x000003E8, 0x4289D68F, 0xFF3900AC, 0x01530000, 0x0000001E, 0x4289D68F, 0xFF3900AC, 0x0153E6A0, 0xFF00001E, 0x4289D68F, 0xFF3900AC, 0x01537C53, 0x00000002, 0x00010229, 0x06890000, 0x0000001E, 0x42366658, 0xFF16007B, 0x002520BA, 0x0000001E, 0x42366658, 0xFF16007B, 0x0025D5E0, 0x000003E8, 0x42366658, 0xFF16007B, 0x00250950, 0x0000001E, 0x42366658, 0xFF16007B, 0x00257C50, 0xFF00001E, 0x42366658, 0xFF17007B, 0x00250000, 0x00000002, 0x0001026F, 0x03500000, 0x0000001E, 0x428D3328, 0xFFCC0011, 0x005B20BA, 0x0000001E, 0x428D3328, 0xFFCC0011, 0x005BD5E0, 0x00000014, 0x428D3328, 0xFFCC0011, 0x005B0950, 0x0000000F, 0x428D3328, 0xFFCC004B, 0x005B7C50, 0x0000000A, 0x42FDFF84, 0xFFFB01F7, 0x00090000, 0x0000001E, 0x42FDFF84, 0xFFFB01F7, 0x0009FFFF, 0x0000001E, 0x42FDFF84, 0xFFFB01F7, 0x00090000, 0x0000001E, 0x42FDFF84, 0xFFFB01F7, 0x0009E6A0, 0xFF00001E, 0x42FDFF84, 0xFFFB01F7, 0x00097C53, 0x00000002, 0x000102B5, 0x043C0000, 0x00000032, 0x42700000, 0x00030006, 0xFFFA20BA, 0x00000028, 0x42700000, 0x00030006, 0xFFFAD5E0, 0x0000001E, 0x424BFFF7, 0x00030006, 0xFFFA0950, 0x00000014, 0x41A4CC7E, 0x00030006, 0xFFFA7C50, 0x00000033, 0x412CCC23, 0x00030006, 0xFFFA0000, 0x00000032, 0x412665BD, 0x00030006, 0xFFFAFFFF, 0x00000032, 0x412665BD, 0x00030006, 0xFFFA0000, 0x00000032, 0x4123328A, 0x00030006, 0xFFFAE6A0, 0xFF000032, 0x412FFF56, 0x00030006, 0xFFFA7C53, 0x00000006, 0x00010301, 0x07BB0000, 0x0000001E, 0x4289332C, 0x00000064, 0x000520BA, 0x0000001E, 0x4289332C, 0x00000065, 0x0006D5E0, 0x0000001E, 0x4289332C, 0x00010063, 0x00290950, 0x0000001E, 0x4289332C, 0x0000002A, 0x00107C50, 0x0000001E, 0x4289332C, 0x0000002A, 0x00100000, 0x000003E8, 0x4289332C, 0x0000002A, 0x0010FFFF, 0x0000001E, 0x4289332C, 0x0000002A, 0x00100000, 0xFF00001E, 0x4289332C, 0x0000002A, 0x0010E6A0, 0xFFFFFFFF, 0x00000000
glabel D_8098875C
 .word func_809856F8
.word func_80985718
.word func_80985738
.word func_80985770
.word func_809857B0
.word func_809857F0
.word func_80985830
.word func_80985C10
.word func_80985C40
.word func_80985C94
.word func_809863BC
.word func_809863DC
.word func_80986430
.word func_80986494
.word func_809864D4
.word func_809868E8
.word func_80986908
.word func_80986948
.word func_80986D40
.word func_80986DC8
.word func_80986E20
.word func_80986E40
.word func_80986EAC
.word func_80986F08
.word func_80986F28
.word func_80986F88
.word func_80986FA8
.word func_80987288
.word func_809872A8
.word func_809872F0
.word func_80987330
glabel D_809887D8
 .word 0x00000000, 0x41200000, 0x00000000
glabel D_809887E4
 .word func_8098764C
.word func_80987658
.word func_80985CE8
glabel Demo_Im_InitVars
 .word 0x00A90400, 0x00000011, 0x00870000, 0x000002FC
.word DemoIm_Init
.word DemoIm_Destroy
.word DemoIm_Update
.word DemoIm_Draw

