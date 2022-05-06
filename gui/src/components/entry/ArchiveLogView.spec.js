import React from 'react'
import 'regenerator-runtime/runtime'
import { render, screen, startAPI, closeAPI, within } from '../conftest.spec'
import ArchiveLogView from './ArchiveLogView'
import EntryContext from './EntryContext'
import userEvent from '@testing-library/user-event'
import { waitFor } from '@testing-library/react'

test('Correctly renders the page', async () => {
  await startAPI('tests.states.uploads.archive_browser_test', 'tests/data/uploads/archive_logs_test', 'test', 'password')
  render(<EntryContext entryId={'1WGSYo1RrGFEIcM17Re4kjHC7k6p'}>
    <ArchiveLogView />
  </EntryContext>)

  // Checking for the page to be successfully loaded and the See more button appears at the bottom
  const seeMoreButton = await screen.findByText(/see more/i)
  expect(seeMoreButton).toBeInTheDocument()
  expect(screen.getByText(/filter logs by level:/i)).toBeInTheDocument()
  expect(screen.getByText(/filter keys by:/i)).toBeInTheDocument()

  // 10 Accordions should be loaded for the given dataset
  const logs = screen.queryAllByTestId('Accordions')
  expect(logs).toHaveLength(10)

  // // Checking for the checkboxes to be checked on the first mount
  const debugBox = screen.getByRole('checkbox', {
    name: /debug/i
  })
  const errorBox = screen.getByRole('checkbox', {
    name: /error/i
  })
  const criticalBox = screen.getByRole('checkbox', {
    name: /critical/i
  })
  const warningBox = screen.getByRole('checkbox', {
    name: /warning/i
  })
  const infoBox = screen.getByRole('checkbox', {
    name: /info/i
  })
  expect(debugBox).toBeChecked()
  expect(errorBox).toBeChecked()
  expect(criticalBox).toBeChecked()
  expect(warningBox).toBeChecked()
  expect(infoBox).toBeChecked()

  // Unchecking the checkbox INFO should re-paint the DOM with only one log and no seeMore button
  userEvent.click(infoBox)
  let updatedLogs = await screen.queryAllByTestId('Accordions')
  expect(updatedLogs).toHaveLength(1)
  expect(seeMoreButton).not.toBeInTheDocument()

  // Re-checking the INFO button should repaint the DOM with the seeMore Button as well as the 10 logs
  userEvent.click(infoBox)
  expect(await screen.findByText(/see more/i)).toBeInTheDocument()
  expect(await screen.queryAllByTestId('Accordions')).toHaveLength(10)

  // Selecting new key from the dropdown menu would add that key to the description of all logs that exist
  await waitFor(async () => expect(await screen.queryByTestId('system_size')).not.toBeInTheDocument())
  const view = await screen.getByTestId('selectOption')
  const butt = await within(view).getByRole('button')
  await userEvent.click(butt)
  await waitFor(async () => expect(await screen.queryByTestId('system_size')).toBeInTheDocument())
  userEvent.click(await screen.queryByTestId('system_size'))
  expect(await screen.getByText(/debug: parsers\/vasp \| undefined/i))

  closeAPI()
})
