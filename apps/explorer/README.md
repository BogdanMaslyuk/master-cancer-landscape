# MCL Explorer frontend

Next.js research interface for Master Cancer Landscape.

## Preferred local run

From the repository root, use the canonical start script:

```powershell
.\scripts\start-frontend.ps1
```

The script sets `NEXT_PUBLIC_MCL_API_URL=http://127.0.0.1:8000`, builds the production frontend automatically when `.next/BUILD_ID` is absent, and then runs `next start`.

The backend should be running separately with:

```powershell
.\scripts\start-backend.ps1
```

Open `http://127.0.0.1:3000`.

## Development mode

For UI work only:

```powershell
cd apps\explorer
$env:NEXT_PUBLIC_MCL_API_URL="http://127.0.0.1:8000"
npm.cmd run dev
```

## Verification

```powershell
npm.cmd run typecheck
npm.cmd run build
```

Or run the full project gate from the repository root:

```powershell
.\scripts\verify.ps1
```
