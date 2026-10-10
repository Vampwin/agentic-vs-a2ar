# สรุปผลทั้งหมดของ "การเทียบ Vina กับ DiffDock ที่ search space เท่ากัน" (ฉบับร่าง 2026-10-09, ยังไม่ได้แก้เอกสารใดๆ)

อ้างอิง: ADR-AIIM-353/354 (Vina blind, DiffDock-L ตัดโปรตีน), ADR-AIIM-355/356 (DiffDock-L กรองด้วย box, DiffDock-Pocket).
**สถานะ:** variant เสริมของ DiffDock-Pocket (rigid, C-alpha centre) ยังรันอยู่ — ช่องที่รอผลมี `[รอ]` กำกับ

## 1. ทำอะไรไปบ้าง

| arm | search space | โปรแกรม | ผ่าน control? | ตีความได้? |
|---|---|---|---|---|
| Vina-box (เดิม, ADR-313) | box 22 Å รอบ site | Vina | 0.40 Å | ใช่ |
| DiffDock-L blind (เดิม, ADR-111/317; 10 samples) | ทั้งโปรตีน | DiffDock-L | 1.92 Å | ใช่ |
| **Vina blind** (A, ใหม่) | box ทั้ง receptor, exhaustiveness 64 | Vina | **0.37 Å** | ใช่ |
| DiffDock-L โปรตีนตัด 12 Å (B, ADR-353) | ใกล้ site | DiffDock-L | **10.1 Å ล้ม** | ไม่ |
| DiffDock-L 100 samples กรองด้วย box (D) | ทั้งโปรตีน แล้วกรอง | DiffDock-L | **5.92 Å ล้ม** (ท่าถูก 10 ท่าอยู่ใน 100 แต่ confidence จัดผิดอันดับ) | ไม่ |
| **DiffDock-Pocket** (C, ติดตั้ง ADR-355) | pocket ที่ centroid ZMA + 6 side chain flexible | DiffDock-Pocket | **0.35 Å** | ใช่ |
| DiffDock-Pocket rigid / C-alpha centre (เสริม) | — | — | [รอ] | [รอ] |

## 2. ผลเทียบแบบ search space เท่ากัน (สาร 15 ตัว, ตรงกัน = RMSD ≤ 2.0 Å)

| คู่ที่เทียบ | ตรงกัน | RMSD มัธยฐาน | หมายเหตุ |
|---|---|---|---|
| Vina-box vs DiffDock-L blind (เดิม ไม่เท่ากัน) | 0/15 | 15.3 Å | ADR-327 |
| **Vina blind vs DiffDock-L blind (blind เท่ากัน)** | **0/15** | **15.9 Å** | ทั้งคู่ผ่าน control |
| **Vina-box vs DiffDock-Pocket (pocket เท่ากัน)** | **0/15** | **12.8 Å** | ใกล้สุด SNP00163 3.4 Å, SNP00006 3.5 Å |
| Vina-box vs DiffDock-L กรอง box | 0/15 | 13.9 Å | ไม่ตีความ (control ล้ม) |

จำนวนที่อยู่ใน pocket (≤ 8 Å จาก centroid ZMA): Vina-box 7, **Vina blind 2**, DiffDock-L blind 11, DiffDock-Pocket 15 (เป็นผลของการให้ pocket), DiffDock-L กรอง box 15 (เป็นผลของการกรอง)

## 3. ข้อสรุปที่ข้อมูลรองรับ

1. **ความไม่ตรงกันของสองโปรแกรมไม่ได้เกิดจาก box**: เมื่อ blind ทั้งคู่ 0/15, เมื่อให้ pocket ทั้งคู่ (DiffDock-Pocket) 0/15
2. **ผลใน pocket ของ Vina ขึ้นกับ box**: เมื่อค้นทั้ง receptor สาร 13/15 ตัวมีท่าดีที่สุดห่าง site 11.9–13.3 Å (ไม่มีท่าใน 10 ท่าที่เก็บอยู่ใน 8 Å ในสาร 13 ตัว); 8 ตัวที่เคยเป็น "off-pocket" อยู่ส่วนนอกของ box (11.9–12.7 Å จาก centroid, box ครึ่งความกว้าง 11 Å); คะแนน blind ดีกว่าคะแนนใน box เกิน 0.5 kcal/mol ใน 8/15 ตัว
3. **สารที่ทำ MD**: SNP00157/SNP00163 ท่าใน site ของ Vina (3.8 และ 1.9 Å จาก centroid) ไม่ใช่ท่าที่ Vina ชอบเมื่อค้นทั้งโปรตีน (ห่าง 11.9 และ 12.0 Å); DiffDock-Pocket วางสองตัวนี้ห่างท่า Vina 9.9 และ 3.4 Å MD ทดสอบท่าใน site ตามที่ ADR-316 ระบุ ไม่ใช่ว่าสารไปอยู่ที่ไหนเมื่อไม่ถูกบังคับ
4. DiffDock-L เจอ site แต่การจัดอันดับด้วย confidence ไม่เสถียร (control: ท่าถูก 10/100 ท่า แต่อันดับ 1 ผิดใน 100 samples; รอบ 10 samples เคยได้ถูก) — rank-1 ของ DiffDock-L รอบเดียวมี noise
5. **ที่ไม่ได้บอก**: ใครถูก, สารจับจริงหรือไม่, อะไรเกี่ยวกับกิจกรรม DiffDock-Pocket ยืดหยุ่น side chain 6 ตัวส่วน Vina แข็ง, รันครั้งเดียว ไม่มี seed, ไม่มี bilayer ทั้งคู่, receptor เป็น 4EIY ที่ตัด BRIL (มี ZMA จับอยู่)

## 4. สิ่งที่ผมแก้/พบที่ต้องแก้ในเอกสาร (เสนอ ยังไม่ได้ทำ)

- **Supplement Section 14 บอก "at 20 samples per complex" ซึ่งผิด**: config YAML ทับ command line จึงได้ 10 samples จริง (ADR-355) ควรแก้เป็น 10
- ประโยค Supplement "Vina searched a 22 Å box ... could not place a ligand elsewhere" ถูกต้องอยู่แล้ว แต่ควรชี้ไปที่ผลที่ใช้ search space เท่ากัน
- §3.8 ย่อหน้า DiffDock: ประโยค "The pose of each survivor is therefore engine-dependent, and the on-pocket result for the shortlisted compounds rests on Vina alone" ต้องปรับ

### ข้อความภาษาอังกฤษที่เสนอ (§3.8, แทรกหลังประโยค "...(Figure 6B-E; Supplementary Table S10)." และแทนประโยคถัดไป)

> Because Vina had been given a 22 Å box and DiffDock the whole protein, the comparison was repeated with matched search spaces (ADR-AIIM-353 to 356). Run over the whole receptor, Vina recovered the ZM241385 crystal pose (0.37 Å) but placed only 2 of the 15 survivors within 8.0 Å of the site; for 13 its best pose lay 11.9 to 13.3 Å from it, so the on-pocket result of the box run depends on the box. Blind against blind, none of the fifteen poses of the two engines agreed within 2.0 Å (median 15.9 Å). DiffDock-Pocket, given the site as the pocket centre and six flexible side chains, reproduced the control to 0.35 Å and placed all fifteen survivors within 1.7 Å of the site centre, yet none of its poses lay within 2.0 Å of the Vina pose in the box (median 12.8 Å; nearest, SNP00163 and SNP00006, 3.4 and 3.5 Å). Two other ways of giving DiffDock-L a pocket failed their own control (a cropped protein, 10.1 Å; 100 blind samples filtered by the box, 5.9 Å) and are not interpreted. The pose of each survivor is therefore engine-dependent even when both engines search the same space; this does not show that the Vina poses are wrong, only that no second method corroborates them.

เพิ่มประโยคเกี่ยวกับ MD ต่อท้ายย่อหน้าที่พูดถึง SNP00157/163 (หรือใน Limits): "Run over the whole receptor, Vina preferred other places for both compounds (11.9 and 12.0 Å from the site), so the simulations tested the in-site poses and say nothing about where the compounds would go unconstrained."

### Figure 6 (ตัวเลือก — ต้องให้คุณเลือก)
- **ก)** คงเดิมแล้วเพิ่ม panel F: จุดต่อสารของ RMSD เทียบ Vina-box ของ DiffDock-L blind (เดิม) กับ DiffDock-Pocket, เส้น 2.0 Å (ใช้ตัวเลขใน `filter_pocket_poses.csv`)
- **ข)** เปลี่ยน panel E (centroid Vina vs DiffDock) เป็น 2×2 ของจำนวนใน pocket
- **ค)** ไม่แก้รูป ใช้ตารางใหม่ใน Supplement แทน (Table S12: ต่อสาร 16 แถว คอลัมน์ centroid/RMSD/confidence ของทุก arm) — เบาที่สุด ผมแนะนำ ค) + legend ระบุว่า panel D/E เป็นการเทียบที่ search space ไม่เท่ากัน

## 5. งานที่เหลือหลังคุณตัดสินใจ

1. รอ variant เสริมของ DiffDock-Pocket (คาด ~22:00 น.) แล้วบันทึกผลใน ADR-AIIM-358
2. แก้บทความ/supplement (tracked changes) ตามที่อนุมัติ, สร้างสำเนาไฮไลต์เหลือง, ตรวจ reject-all
3. เพิ่ม `docking_matched_search_space_2026-10-09` และ `docking_matched_filter_pocket_2026-10-09` ใน EXTENSION_TREE ของ `build_deposit.py` (**ต้องอนุมัติ**; ไม่เอาน้ำหนักโมเดลและ `.cache`), build/amend commit, อัปเดต hash/จำนวนไฟล์ในบทความ
4. รัน G8 ใหม่ + new_number_audit ตัวเลขใหม่ (15.9, 12.8, 0.35, 0.37, 1.7, 13/15 ฯลฯ) — G8 เป็นตัวจับ literal ตัวเลขใหม่จะไม่ถูกตรวจจนกว่าจะเพิ่มตรวจเอง
5. ใส่ README ที่บอกว่าสคริปต์ใช้พาธ absolute และไม่ได้ฝาก DiffDock-Pocket (ติดตั้งตาม ADR-355)

## 6. ที่ต้องให้คุณตัดสินใจ

- อนุมัติข้อความด้านบน (หรือขอแก้)
- เลือกตัวเลือกรูป ก/ข/ค
- อนุมัติแก้ `build_deposit.py` เพื่อใส่สองโฟลเดอร์ใหม่
- ยังต้อง `git push` เองและตั้ง repository เป็น public
