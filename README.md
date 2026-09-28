# jouzoukettei — 有機構造決定パズル

京大化学の有機（大問3の構造決定）で使う「条件から候補を絞り込む推論」を反復訓練するブラウザゲーム。
仕様は Claude Docs の「有機構造決定パズル 仕様書」。いまは MVP（絞り込み型・C₅以下の脂肪族・20問）。

## 遊び方

`dist/index.html` をブラウザで開くだけ。サーバーもネット接続もいらない（フォントだけ Google Fonts から読む）。
`dist/` は GitHub Actions が自動で生成してコミットする。

1. 問題を選ぶ（難易度 ★1〜★3、または「おまかせ」）
2. 手がかりカードをタップして結果を見る。1枚ごとに −10 点
3. 矛盾する候補を自分で消す（左右スワイプか「消す」）。自動では消えない
4. 残り1つで「確定」。誤答なら矛盾したカード、根拠なしで消した候補は消すべきだったカードが出る
5. 間違えたカードの種類が「復習リスト」に溜まり、おまかせ出題で多めに出る

23時〜5時はプレイできない（仕様書の運用ルール）。

## 構成

| パス | 役割 |
| --- | --- |
| `problems/narrow.json` | 問題データ（仕様書の JSON 形式＋`level`） |
| `problems/names.json` | 結果画面に出す物質名 |
| `src/chem.js` | 手がかりカードの判定ロジック（SMARTS・酸化/加水分解/オゾン分解・不斉炭素・シス-トランス） |
| `scripts/test-rules.js` | 判定ロジックの単体テスト（教科書の典型例） |
| `scripts/validate.js` | 問題の自動検証 |
| `scripts/build.js` | 検証 → 構造式SVGと判定表を前計算 → `dist/index.html` を出力 |
| `src/index.html` | ゲーム画面のテンプレート |

化学処理（RDKit）はビルド時に Node で行い、ブラウザには前計算した構造式SVGと判定結果だけを埋め込む。
仕様書ではブラウザで RDKit.js を CDN から読む案だったが、wasm 読み込みが不要になり、スマホでも一瞬で開けてオフラインでも動くのでこちらにした。

## 問題を足す

```sh
npm install
npm test          # 判定ロジックのテスト + 全問題の検証
npm run build     # dist/index.html を作り直す
```

検証スクリプトが確認すること（1つでも落ちた問題は収録しない）:

1. 全候補の分子式が `formula` と一致する
2. 候補に重複がない（canonical SMILES で比較）
3. 全手がかりを適用すると、候補が `answer` の1つだけ残る
4. どの手がかりも、単独で候補を1つ以上消す
5. 手がかりの `result` が、正解の構造に実際にカードを当てた結果と一致する
6. カードが適用できない構造（例: C=C を持つ候補に KMnO₄）を含まない

## 手がかりカード

| card | 画面の名前 | result |
| --- | --- | --- |
| `silver_mirror` | 銀鏡反応 | true / false |
| `fehling` | フェーリング液 | true / false |
| `iodoform` | ヨードホルム反応 | true / false |
| `sodium` | 金属Na | true / false |
| `nahco3` | NaHCO₃ | true / false |
| `fecl3` | FeCl₃ | true / false |
| `bromine` | 臭素水 | true / false |
| `kmno4` | KMnO₄ 酸化 | 生成物 SMILES の配列（`[]` = 酸化されない） |
| `hydrolysis` | 加水分解 | 生成物 SMILES の配列（`[]` = 加水分解されない） |
| `ozonolysis` | オゾン分解 | 生成物 SMILES の配列（`[]` = 反応しない） |
| `chiral` | 不斉炭素 | 個数 |
| `cis_trans` | シス-トランス異性 | true / false |

ゲーム画面にはカードの結果しか出さない。「その結果から何が確定するか」は仕様書の手がかり表の右端の列に自分で書く。
