export default function ReloadSnack() {
  return <Snackbar
    anchorOrigin={{
      vertical: 'bottom',
      horizontal: 'left',
    }}
    open
    message={<span id="message-id">There is a new version. Please press your browsers Reload (or even Shift+Reload) button.</span>}
/>
}